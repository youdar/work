from __future__ import division
from iotbx.pdb.multimer_reconstruction import multimer
import mmtbx.refinement.minimization_ncs_constraints
import mmtbx.monomer_library.pdb_interpretation
from iotbx.pdb import format_MTRIX_pdb_string
import mmtbx.ncs.ncs_utils as nu
import mmtbx.monomer_library.server
from libtbx import adopt_init_args
from libtbx.utils import null_out
from iotbx.pdb import fetch
from iotbx import mtz
import mmtbx.f_model
import mmtbx.utils
import iotbx.ncs
import iotbx.cif
import iotbx.pdb
import tempfile
import getpass
import time
import os
import sys

__author__ = 'Youval Dar'

# Parameters needed for ADP energies
adp_restraints_master_params = iotbx.phil.parse("""\
  iso {
    use_u_local_only = False
      .type = bool
    sphere_radius = 5.0
      .type = float
    distance_power = 1.69
      .type = float
    average_power = 1.03
      .type = float
    wilson_b_weight_auto = False
      .type = bool
    wilson_b_weight = None
      .type = float
    plain_pairs_radius = 5.0
      .type = float
    refine_ap_and_dp = False
      .type = bool
  }
""")

class ncs_refine_test(object):
  def __init__(self,
               working_paths,
               use_geometry_restraints    = True,
               real_data                  = True,
               r_work_reported_pdb_ncs    = -1,
               r_free_reported_pdb_ncs    = -1,
               r_work_calc_pdb_ncs        = -1,
               r_free_calc_pdb_ncs        = -1,
               initial_r_work             = -1,
               initial_r_free             = -1,
               r_work_final               = -1,
               r_free_final               = -1,
               resolution                 = -1,
               pdb_header_year            = -1,
               pdb_code                   = '____',
               use_strict_ncs             = None,
               print_during_refinement    = False,
               sf_and_grads_algorithm     = 'fft',
               target_name                = 'ml',
               use_hd                     = False,
               rotation_matrices          = [],
               translation_vectors        = []):
    """
    Run strict NCS refinment tests
    1) Refine without using strict NCS
    2) The same refinement using strict NCS

    Argument:
    use_geometry_restraints: (bool) Use geometry restraints in refinement
    real_data: (bool) False for testing, True for real data

    Usage examples:
    >>>python test_ncs_refinement.py 1vcr
    >>>python test_ncs_refinement.py 1vcr -quiet

    Output example:
#  PDB code |Reported in PDB  | Calc from NCS   |  ASU initial    |   ASU final     | Res.   |   NCS   |  Solvent |    Data      | Year   | Use  | Geo. | Trans. | time (sec)
#           | r-work | r-free | r-work | r-free | r-work | r-free | r-work | r-free |        | copies  | fraction | completeness |        | NCS  | Rest | refine |
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
     1vcr   | 0.3790 | 0.3530 | 0.4292 | 0.4486 | 0.4242 | 0.4382 | 0.3731 | 0.4521 |  9.50  |    5    |   0.94   |     1.00     |  2004  | True | True |  True  |  236
    """
    # Start time
    self.start_time = time.time()
    # creates self.xxx for all arguments
    adopt_init_args(self, locals())
    self.pdb_asu_dir = working_paths.pdb_asu_dir
    self.pdb_dir = working_paths.pdb_dir
    self.mtz_dir = working_paths.mtz_dir

  def get_structure_factors(self):
    """
    Get f_obs and r_free_flags
    From cif file if available
    """
    self.f_obs = None
    self.r_free_flags = None
    fobs = ["FOBS,SIGFOBS",'FOBS','FOBS,PHIM',
            "F(+),SIGF(+),F(-),SIGF(-)","F(+),F(-)"]
    iobs = ["IOBS,SIGIOBS",'IOBS','IOBS,PHIM',
            'I(+),SIGI(+),I(-),SIGI(-)']
    fn = self.pdb_code + '.mtz'
    fn = os.path.join(self.mtz_dir,fn)
    mtz_object = mtz.object(file_name=fn)
    miller_arrays = mtz_object.as_miller_arrays()
    # print miller_arrays[0].completeness()
    for ma in miller_arrays:
      ls = ma.info().label_string()
      if (ls in fobs) or ma.is_xray_amplitude_array():
        # Consider using Bijvoet mates
        ma = ma.average_bijvoet_mates()
        self.f_obs = abs(ma)
      elif ls == "R-free-flags":
        self.r_free_flags = abs(ma)
      elif not self.f_obs and ((ls in iobs) or ma.is_xray_intensity_array()):
        # Consider using Bijvoet mates
        ma = ma.average_bijvoet_mates()
        # convert i_obs to f_obs
        self.i_obs = ma
      elif not self.r_free_flags and ls == "R-free-flags-1":
        self.r_free_flags = abs(ma.french_wilson(log=null_out()))
    if not self.f_obs and self.i_obs:
      self.f_obs = abs(self.i_obs.french_wilson(log=null_out()))

    if self.f_obs:
      # self.f_obs.show_summary()
      if self.r_free_flags:
        self.f_obs, self.r_free_flags = self.f_obs.common_sets(
          self.r_free_flags)
        self.r_free_flags = self.make_r_free_boolean(self.r_free_flags)
      else:
        self.r_free_flags = self.f_obs.generate_r_free_flags()
      # Data completeness: Fraction of unmeasured reflections within the
      # [d_min, d_max] range,where d_min and d_max are highest and lowest
      # resolution of data set correspondingly.
      self.completeness = self.f_obs.array().completeness()
    else:
      raise RuntimeError("Missing amplitude array.")

  def make_r_free_boolean(self,r_free_flags):
    '''
    Convert r_free_flag from any of the conventions possible
    to a boolean

    Posiblle convention for free and working set flags:
    CCP4     assigns the flag r_free_flags to be 0 for the free set and 1,
             ...n-1 for the working set.
    XPLOR    assigns the flag TEST to be 1 for the free set and 0 for the
             working set.
    CNS      assigns the flag TEST to be 1 for the free set and 0,2,...n-1 for
             the working set.
    SHELX    assigns a flag with -1 for the free set and 1 for the working set.
    TNT      assigns a flag with 0 to indicate the free set.
    '''
    flag_value = iotbx.reflection_file_utils.guess_r_free_flag_value(
      miller_array =r_free_flags)
    if flag_value is None: return None
    else:
      return r_free_flags.customized_copy(data=r_free_flags.data()==flag_value)

  def get_pdb_published_data(self,pdb_inp, xrs_pdb):
    """
    Get the r-work value reported in PDB
    Calculate r-work from a single NCS copy (from f-model)

    Arguments:
    pdb_inp: PDB input object of the single NCS
    xrs_pdb: x-ray structure object for a single NCS
    """
    self.pdb_header_year = pdb_inp.extract_header_year()
    pio = pdb_inp.get_r_rfree_sigma()
    self.r_work_reported_pdb_ncs = pio.r_work
    self.r_free_reported_pdb_ncs = pio.r_free
    self.resolution = pio.resolution
    # calculate r-work, r-free from pdb
    _,self.r_work_calc_pdb_ncs,self.r_free_calc_pdb_ncs = \
      self.get_f_model(xrs_pdb)

  def process_pdb_files(self):
    self.get_structure_factors()
    # Collect r-work for the single NCS copy
    # self.get_pdb_published_data(pdb_inp=pdb_inp,xrs_pdb=xrs_ncs)


  def get_pdb_published_data(self):
    pass

  def refinement_loop(self,
                      n_macro_cycle = 10,
                      sites = True,
                      u_iso = True,
                      transformations = True,
                      finite_grad_differences_test = False,
                      r_work_target = 0.01,
                      use_strict_ncs = True,
                      write_refined_pdb = False):
    """
    Refine the complete ASU, using strict-NCS refinement

    If more than one refinement method (sites, u_iso, transformations) is
    True, refinement process will alternate between those methods

    Args:
      n_macro_cycle: (int) number of refinement cycles
      sites: (bool) refinement using sites
      u_iso: (bool) refinement using b-factors (ADP)
      finite_grad_differences_test: (bool) Compare calculated gradient to
        estimated one
      r_work_target: (float) stop refinement when r_work <= r_work_target
      use_strict_ncs: (bool) When True, use strict NCS refinement
      write_refined_pdb: (bool) When True, write refined pdb file in current
        directory

    """
    pass

class Set_paths(object):


  def __init__(self):
    # check if working at LBL local machine
    osType = sys.platform
    self.pdb_dir = ''
    if osType.startswith('win'):
      self.work_env = 'win'
      self.pdb_asu_dir = r'C:\Phenix\Dev\Work\work\NCS\data\pdb_asu'
      self.pdb_dir = r'C:\Phenix\Dev\Work\work\NCS\data\pdb'
      self.cif_dir = r'C:\Phenix\Dev\Work\work\NCS\data\cif'
      self.mtz_dir = r'C:\Phenix\Dev\Work\work\NCS\data\mtz'
      self.file_index =r'C:\Phenix\Dev\Work\work\NCS\data\pdb_157_file_list.txt'
    else:
      try:
        root = '/net/cci/youval/Work/work/NCS/data/'
        self.pdb_asu_dir = root + 'pdb_asu'
        self.pdb_dir = root + 'pdb'
        self.cif_dir = root + 'cif'
        self.mtz_dir = root + 'mtz'
        self.file_index = root + 'pdb_157_file_list.txt'
        self.work_env = 'lbl'
      except KeyError:
        self.work_env = 'other'
    # parameters for ADP energies
    adp_restraints_params  = adp_restraints_master_params.extract()
    self.iso_restraints = adp_restraints_params.iso
    self.pdb_codes = open(self.file_index,'r').read().splitlines()

def run():
  print 'start'
  print '='*50
  paths = Set_paths()
  # list test to perform
  test_param = [[False,True,False],
                [True,True,False],
                [True,True,True]]
  processed_files = set()
  for pdb_code in paths.pdb_codes:
    t1 = os.path.isfile(os.path.join(paths.pdb_dir,pdb_code + '.pdb'))
    t2 = os.path.isfile(os.path.join(paths.pdb_asu_dir,pdb_code + '-asu.pdb'))
    t3 = os.path.isfile(os.path.join(paths.mtz_dir,pdb_code + '.mtz'))
    have_all_files = t1 and t2 and t3
    if not have_all_files: continue
    processed_files.add(pdb_code)
    for i, (p1,p2,p3) in enumerate(test_param):
      run_test(pdb_code=pdb_code,
               use_strict_ncs=p1,
               use_geometry_restraints = p2,
               refine_transformations = p3,
               working_paths = paths)
  print '='*50
  print 'Number of files processed: ',len(processed_files)
  not_processed_files =  set(wp.pdb_codes) - processed_files
  print 'Number of file not processed:',len(not_processed_files)

def run_test(pdb_code,
             use_strict_ncs,
             use_geometry_restraints,
             refine_transformations,
             working_paths):

  # Run tests without NCS
  test_obj = ncs_refine_test(
    use_geometry_restraints=use_geometry_restraints,
    real_data=True,
    pdb_code=pdb_code,
    working_paths=working_paths)

  test_obj.process_pdb_files()

  # Run refinement
  test_obj.refinement_loop(
    n_macro_cycle=5,
    r_work_target=0.1,
    sites=True,
    u_iso=True,
    transformations=refine_transformations,
    use_strict_ncs=use_strict_ncs)


if __name__=='__main__':
  run()
