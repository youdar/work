from __future__ import division
import tempfile
import getpass
import time
import os
import sys

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
import iotbx.cif
import iotbx.pdb


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
               print_during_refinement    = True,
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
    # check if working at LBL local machine
    osType = sys.platform
    self.pdb_dir = ''
    if osType.startswith('win'):
      self.work_env = 'win'
    else:
      try:
        self.pdb_dir = os.environ["PDB_MIRROR_PDB"]
        self.cif_dir = os.environ["PDB_MIRROR_STRUCTURE_FACTORS"]
        self.mtz_dir = '/net/cci-filer2/raid1/share/pdbmtz/mtz_files'
        self.work_env = 'lbl'
      except KeyError:
        self.work_env = 'other'
    # parameters for ADP energies
    adp_restraints_params  = adp_restraints_master_params.extract()
    self.iso_restraints = adp_restraints_params.iso


  def process_pdb_and_cif_files(self, args):
    """
    Fetch pdb file when not working at LBL, if does not exist in working
    directory

    Argument:
    args: (str) 4 letters pdb code or a list where the first variable
          is (str) 4 letters pdb code
    """
    if type(args) == list:
      file_name = args[0]
    else:
      file_name = args
    assert len(file_name) == 4
    assert type(file_name) == str
    # Read pdb file when not working at LBL
    self.pdb_code = file_name.lower()
    self.pdb_file_name = self.pdb_code + '.pdb'
    self.cif_file_name = self.pdb_code + '-sf.cif'
    # Get pdb file
    if os.path.isfile(self.pdb_file_name):
      self.full_path_pdb = os.path.realpath(self.pdb_file_name)
      print 'Using pdb file from local machine (not from MIRROR)'
    else:
      self.full_path_pdb = self.get_full_path(data_type='pdb')
    # Get cif file
    if os.path.isfile(self.cif_file_name):
      self.full_path_cif = os.path.realpath(self.cif_file_name)
      print 'Using cif file from local machine (not from MIRROR)'
    else:
      self.full_path_cif = self.get_full_path(data_type='xray')
    # Collect Fobs and r_free_flags from cif file
    self.get_structure_factors()
    # Process PDB file
    self.process_pdb_file()

  def process_pdb_file(self):
    """
    Process file
    - Check if it contains a single NCS copy
    - Build xray structure, PDB Hierarchy
    - Get Crystal symmetry
    - Collect MTRIX transformation info
    - Collect NCS atoms selection
    - Get geometry restraints manager (grm)
    - Get reported r-work and r-work of a single NCS

    Arguments:
    build_unit_cell: (bool) leave original crystal symmetry or
                     create a unit cell
    """
    m = multimer(
      file_name=self.full_path_pdb,
      round_coordinates=False,
      reconstruction_type='cau',
      error_handle=True,eps=1e-2)
    #
    crystal_symmetry = self.f_obs.crystal_symmetry()
    pdb_inp = iotbx.pdb.input(file_name=self.full_path_pdb)
    ph_ncs = pdb_inp.construct_hierarchy()
    xrs_ncs = pdb_inp.xray_structure_simple()
    crystal_sym_test = crystal_symmetry.\
      is_similar_symmetry(xrs_ncs.crystal_symmetry())
    if not crystal_sym_test:
      print 'Warning :crystal_symmetry of PDB is different than the one from' \
            ' Fobs. Using crystal_symmetry from Fobs.'
      xrs_ncs = pdb_inp.xray_structure_simple(crystal_symmetry=crystal_symmetry)
    pdb_inp_asu = m.assembled_multimer.as_pdb_input()
    xrs_asu = pdb_inp_asu.xray_structure_simple(
      crystal_symmetry = crystal_symmetry)
    # force ASU none-rounded coordinates into xray structure
    xrs_asu.set_sites_cart(m.sites_cart())
    assert crystal_symmetry.is_similar_symmetry(xrs_asu.crystal_symmetry())
    # Collect r-work for the single NCS copy
    self.get_pdb_published_data(pdb_inp=pdb_inp,xrs_pdb=xrs_ncs)
    # Set attributes
    self.ncs_hierarchy = ph_ncs
    self.crystal_symmetry = crystal_symmetry

    # print pdb_asu to data folder
    path_pdb_asu='/net/cci-filer2/raid1/home/youval/Work/work/NCS/data/pdb_asu'
    fn_dest = os.path.join(path_pdb_asu,self.pdb_code + '-asu.pdb')
    m.write(pdb_output_file_name=fn_dest,
            crystal_symmetry=self.crystal_symmetry)

    self.rotation_matrices = m.rotation_matrices
    self.translation_vectors = m.translation_vectors
    # keep transformation data for comparison at the end of refinement
    self.init_rotation_matrices = m.rotation_matrices
    self.init_translation_vectors = m.translation_vectors
    selection_str = 'chain ' + ' or chain '.join(m.ncs_unique_chains_ids)
    self.ncs_atom_selection = m.assembled_multimer.\
      atom_selection_cache().selection(selection_str)
    if self.ncs_atom_selection.count(True) == 0:
      print 'Atom selection issues'
    # Get geometry restraints manager (grm)
    self.grm = None
    if self.print_during_refinement:
      print 'Time till start GRM calc: {}'.format(self.time_from_start())
    if self.use_geometry_restraints:
      pdb_str = m.assembled_multimer.as_pdb_string(
        crystal_symmetry=crystal_symmetry)
      self.grm = nu.get_restraints_manager(pdb_string=pdb_str)

    if self.print_during_refinement:
      print 'Time till end GRM calc: {}'.format(self.time_from_start())
    # Get f_model and initial r_work
    self.fmodel, self.initial_r_work, self.initial_r_free = \
      self.get_f_model(xrs=xrs_asu)
    if self.print_during_refinement:
      print 'Time till end of Fmodel calc: {}'.format(self.time_from_start())
    #
    self.get_solvent_fraction(xrs_asu)

  def get_solvent_fraction(self, xrs_asu):
    """
    The fraction of the weight, of the unit cell, contributed
    by solvent
    """
    import mmtbx.masks
    self.solvent_fraction = mmtbx.masks.asu_mask(
      xray_structure=xrs_asu,
      d_min=1).asu_mask.contact_surface_fraction

  def get_pdb_published_data(self,pdb_inp, xrs_pdb):
    """
    Get the r-work value reported in PDB
    Calculate r-work from a single NCS copy (from f-model)

    Arguments:
    pdb_inp: PDB input object of the single NCS
    xrs_pdb: x-ray structure object for a single NCE
    """
    self.pdb_header_year = pdb_inp.extract_header_year()
    pio = pdb_inp.get_r_rfree_sigma()
    self.r_work_reported_pdb_ncs = pio.r_work
    self.r_free_reported_pdb_ncs = pio.r_free
    self.resolution = pio.resolution
    # calculate r-work, r-free from pdb
    _,self.r_work_calc_pdb_ncs,self.r_free_calc_pdb_ncs = \
      self.get_f_model(xrs_pdb)

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
    if self.full_path_cif:
      # miller_arrays = iotbx.cif.reader(
      #   file_path=self.full_path_cif).\
      #   as_miller_arrays(force_symmetry=True)
      miller_arrays = self.get_miller_arrays()
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
          self.f_obs = abs(ma.french_wilson(log=null_out()))
        elif not self.r_free_flags and ls == "R-free-flags-1":
          self.r_free_flags = abs(ma.french_wilson(log=null_out()))
    else:
      raise RuntimeError("No cif file.")

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

  def get_miller_arrays(self):
    """
    convert cif of the format 'r_name_sf.ent.gz' to mtz file
    creates mtz file with crystal symmetry in current folder
    """
    have_mtz_file = self.full_path_cif[-3:] == 'mtz'
    if not have_mtz_file:
      from libtbx import easy_run
      tempmtzFile = tempfile.NamedTemporaryFile().name + '.mtz'
      cmd_list = []
      cmd_list.append('phenix.cif_as_mtz')
      cmd_list.append(self.full_path_cif)
      cmd_list.append('--output-file-name={}'.format(tempmtzFile))
      cmd_list.append('--symmetry={}'.format(self.full_path_pdb))
      cmd_list.append("--merge")
      cmd_list.append("--remove-systematic-absences")
      cmd_list.append("--map-to-asu")
      cmd_list.append("--ignore-bad-sigmas")
      cmd_list.append("--extend-flags")
      cmd = ' '.join(cmd_list)
      r = easy_run.go(cmd)
      # NOTE !!! in windows r does not returns the errors as expected
      tmp = [x for x in r.stdout_lines if '--' in x]
      for x in tmp:
        tmp2 = [xx for xx in x.split() if '--' in xx]
        s  = '--incompatible_flags_to_work_set'
        if s in tmp2 and s not in cmd:
          cmd_list.append(s)
          cmd = ' '.join(cmd_list)
          easy_run.go(cmd)
    else:
      tempmtzFile = self.full_path_cif

    try:
      mtz_object = mtz.object(file_name=tempmtzFile)
      # Process the mtz_object
      miller_arrays = mtz_object.as_miller_arrays()
    except:
      miller_arrays = None

    if not have_mtz_file:
      # cleanup - delete temp.mtz
      os.remove(tempmtzFile)
    return miller_arrays


  def get_f_model(self, xrs):
    """
    Argument:
    xrs: x-ray structure factor object

    Return: fmodel,r_work
    fmodel: f-model object
    r_work: (float) r-work derived from the f-model
    """
    params = mmtbx.f_model.sf_and_grads_accuracy_master_params.extract()
    params.algorithm = self.sf_and_grads_algorithm
    fmodel = mmtbx.f_model.manager(
      f_obs                        = self.f_obs,
      r_free_flags                 = self.r_free_flags,
      xray_structure               = xrs,
      sf_and_grads_accuracy_params = params,
      target_name                  = self.target_name)
    # Filter by resolution
    # fmodel = fmodel.resolution_filter(d_min = res)
    # Update model for Bulk-solvent and overall scaling
    # when working with real data
    if self.real_data:
      fmodel.update_all_scales(update_f_part1_for=None)
    r_work = fmodel.r_work()
    r_free = fmodel.r_free()
    return fmodel.deep_copy(), r_work, r_free

  def refinement_loop(self,
                      n_macro_cycle = 10,
                      sites = True,
                      u_iso = True,
                      transformations = True,
                      finite_grad_differences_test = False,
                      r_work_target = 0.01,
                      use_strict_ncs = True,
                      print_during_refinement=False):
    """
    Refine the complete ASU, using strict-NCS refinement

    If more than one refinement method (sites, u_iso, transformations) is
    True, refinement process will alternate between those methods

    Arguments:
    n_macro_cycle: (int) number of refinement cycles
    sites: (bool) refinement using sites
    u_iso: (bool) refinement using b-factors (ADP)
    finite_grad_differences_test: (bool) Compare calculated gradient to
                                  estimated one
    r_work_target: (float) stop refinement when r_work <= r_work_target
    use_strict_ncs: (bool) When True, use strict NCS refinement

    """
    self.use_strict_ncs = use_strict_ncs
    # print only if the not in quiet mode
    print_during_refinement = \
      print_during_refinement and self.print_during_refinement
    if n_macro_cycle > 100:
      raise IOError('To many refinement cycles')
    if self.print_during_refinement:
      print 'Time till refinement cycle start: {}'.format(self.time_from_start())
    # prepare for alternate refinement methods
    refinement_method = self.use_refinement_method(sites,u_iso,transformations)
    n_macro_cycle *= [sites,u_iso,transformations].count(True)
    self.refinement_time_start = time.time()
    for macro_cycle in xrange(n_macro_cycle):
      refinement_method = self.iterate_refine_method(
        refinement_method,macro_cycle)
      self.sites = refinement_method['sites'][2]
      self.u_iso = refinement_method['u_iso'][2]
      self.transformations = refinement_method['transformations'][2]
      if transformations and not self.rotation_matrices: continue
      data_weight = None
      if(self.use_geometry_restraints):
        data_weight = nu.get_weight(self)
      minimized = mmtbx.refinement.minimization_ncs_constraints.lbfgs(
        fmodel                       = self.fmodel,
        rotation_matrices            = self.rotation_matrices,
        translation_vectors          = self.translation_vectors,
        ncs_atom_selection           = self.ncs_atom_selection,
        finite_grad_differences_test = finite_grad_differences_test,
        geometry_restraints_manager  = self.grm,
        data_weight                  = data_weight,
        max_iterations               = 35,
        refine_sites                 = self.sites,
        refine_u_iso                 = self.u_iso,
        refine_transformations       = self.transformations,
        use_strict_ncs               = self.use_strict_ncs,
        iso_restraints               = self.iso_restraints,
        use_hd                       = self.use_hd)
      # For real data: Update model for Bulk-solvent and overall scaling
      if self.real_data:
        self.fmodel.update_all_scales(update_f_part1_for=None)
      if transformations:
        self.rotation_matrices = minimized.rotation_matrices
        self.translation_vectors = minimized.translation_vectors
      if print_during_refinement:
        refine_type = 'adp'*self.u_iso + 'sites'*self.sites \
                      + 'transforms'*self.transformations
        outstr = "  macro_cycle {0:3} ({1:10})   r_factor: {2:6.4f}   " + \
              finite_grad_differences_test * \
              "finite_grad_difference_val:{3:.4f}"
        print outstr.format(
          macro_cycle, refine_type,self.fmodel.r_work(),
          minimized.finite_grad_difference_val)
      if self.fmodel.r_work()<=r_work_target: break
    self.r_work_final = self.fmodel.r_work()
    self.r_free_final = self.fmodel.r_free()

    # save results to data folders
    dest = None
    if not self.use_strict_ncs:
      dest = '/net/cci-filer2/raid1/home/youval/Work/work/NCS/data/refinement_no_ncs/'
    elif transformations:
      dest = '/net/cci-filer2/raid1/home/youval/Work/work/NCS/data/refinement_ncs_transform/'
    else:
      dest = '/net/cci-filer2/raid1/home/youval/Work/work/NCS/data/refinement_ncs/'
    self.write_refined_pdb_file(
      destination=dest,
      as_ncs = self.use_strict_ncs,
      transforms = transformations)


  def write_refined_pdb_file(self,
                             destination = '',
                             as_ncs = True,
                             transforms = False):
    """
    Writes the refined NCS or ASU in path and file name specified by 'destination'

    Argument:
    dest: (str) Path and file name for output
    as_ncs: (bool) when True, save as NCS. When False, save as ASU
    """
    if destination[-4:] != '.pdb':
      s  = as_ncs * '-ncs' + (not as_ncs)*'-asu' + transforms*'-trnsfrm'
      destination += self.pdb_code + s + '-refined.pdb'
    temp_xrs = self.fmodel.xray_structure.select(self.ncs_atom_selection)
    xyz = temp_xrs.sites_cart()
    self.ncs_hierarchy.atoms().set_xyz(xyz)
    mtrix_str = format_MTRIX_pdb_string(
      self.rotation_matrices,self.translation_vectors)

    of = open(destination, "w")
    outstr1 = [
      x for x in temp_xrs.xray_structure().as_pdb_file().splitlines()
      if x.startswith('CRYS') or x.startswith('SCALE')]
    outstr1 = '\n'.join(outstr1)
    print >> of, outstr1
    print >> of, mtrix_str
    print >> of, self.ncs_hierarchy.as_pdb_string()
    of.close()


  def use_refinement_method(self,sites,u_iso,transformations):
    """
    create a dictionary of refinement methods to alternate over

    refine_method[method] = [use (bool), turn (int), initial value (bool)
    Each method that should be used will be set to True on its turn, one at
    the time
    """
    if [sites, u_iso, transformations].count(True) == 0:
      raise IOError('Refinement method error')
    types = ('sites','u_iso','transformations')
    use = (sites,u_iso,transformations)
    i = 0
    method = {}
    for (k,v) in zip(types,use):
      if v:
        method[k] = [True,i,False]
        i += 1
      else:
        method[k] = [False,-1,False]
    return method

  def iterate_refine_method(self,refine_dict,cycle_num):
    """
    Change one of the refinement methods to True, according to the cycle_num
    and the methods that should be used
    """
    use_methods = [k for (k,v) in refine_dict.iteritems() if v[0]]
    n = cycle_num % len(use_methods)
    for k in use_methods:
      if refine_dict[k][1] == n:
        refine_dict[k][2] = True
      else:
        refine_dict[k][2] = False
    return refine_dict

  def cleanup(self):
    """
    Cleanup temporary files and folders created during test
    """
    files_to_delete = ['data_{}.mtz'.format(self.pdb_code)]
    for fn in files_to_delete:
      if os.path.isfile(fn): os.remove(fn)
    os.chdir(self.current_folder)

  def get_full_path(self, data_type='pdb'):
    '''(str) -> str
    returns the full path a pdb file name in LBL PDB_MIRROR_PDB

    Argument:
    file_type: (str) the type of file we are looking for
               'pdb' pdb or 'xray' for cif

    Returns:
    full_path: (str) full path of file in LBL PDB_MIRROR_PDB
    '''
    full_path = None
    if self.work_env == 'lbl':
      # Use file directly from LBL MIRRORs
      if data_type=='pdb':
        mirror_dir = self.pdb_dir
        s=-12; e=-8  # start and end of 4 letter pdb code in record
      elif data_type=='xray':
        mirror_dir = self.cif_dir
        s=-14; e=-10  # start and end of 4 letter pdb code in record
        full_path = os.path.join(self.mtz_dir,self.pdb_code + '.mtz')
      else: raise IOError('Wrong file format')
      if (not full_path) or (not os.path.isfile(full_path)):
        pdb_files = open(os.path.join(mirror_dir, "INDEX"), "r").readlines()
        # assuming that the 4 letter PDB code is in [s:e] of the records
        files = [x.strip() for x in pdb_files if x[s:e]==self.pdb_code]
        if len(files)!=1:
          outstr = 'Found {0} number of {1} in PDB_MIRROR_PDB index'
          print outstr.format(len(files),self.pdb_code)
        else:
          full_path = os.path.join(mirror_dir,files[0])
    else:
      # Fetch file to local current folder
      fetched_file = fetch.get_pdb (
        id=self.pdb_code, data_type=data_type,
        mirror='rcsb',quiet=True,log=sys.stdout)
      full_path = os.path.realpath(fetched_file)
    return full_path

  def set_working_path(self,path=None):
    """
    Change working directory to avoid littering
    """
    self.current_folder = os.getcwd()
    username = getpass.getuser()
    if username.lower() == 'youval':
      if path:
        self.tempdir = path
      else:
        osType = sys.platform
        if osType.startswith('win'):
          self.tempdir = (r'C:\Phenix\Dev\Work\work\NCS\junk\pdb_test')
        else:
          self.tempdir = ('/net/cci/youval/Work/work/NCS/junk/pdb_test')
      os.chdir(self.tempdir)

  def __repr__(self):
    """
    Print refinement and PDB info summery
    """
    # File and R-factors summery
    data_in = [self.pdb_code,
               self.r_work_reported_pdb_ncs,
               self.r_free_reported_pdb_ncs,
               self.r_work_calc_pdb_ncs,
               self.r_free_calc_pdb_ncs,
               self.initial_r_work,
               self.initial_r_free,
               self.r_work_final,
               self.r_free_final,
               self.resolution,
               len(self.rotation_matrices) + 1,
               self.solvent_fraction,
               self.completeness,
               self.pdb_header_year,
               str(self.use_strict_ncs),
               str(self.use_geometry_restraints),
               str(self.transformations),
               self.time_from_start()]
    # data
    data_out = []
    # make sure we do not have non numerical r-values
    for i,x in enumerate(data_in):
      if (not x) and i>0 and i<12:
        data_out.append(-1)
      else: data_out.append(x)
    sep = '\n#' + '-'*172 + '\n'
    s = '  {0:^10}|{1:^8.4f}|{2:^8.4f}|{3:^8.4f}|{4:^8.4f}|{5:^8.4f}'
    s += '|{6:^8.4f}|{7:^8.4f}|{8:^8.4f}|{9:^8.2f}|{10:^9}|{11:^10.2f}'
    s += '|{12:^14.2f}|{13:^8}|{14:^6}|{15:^6}|{16:^8}|{17:^8}'
    data_str = s.format(*data_out)
    # first title row
    t1 = '# {0:^10}|{1:^16} |{2:^16} |{3:^16} |{4:^16} |{5:^7} | {6:^7}'
    t1 += ' | {7:^9}|{8:^13} |{9:^7} |{10:^6}|{11:^6}|{12:^8}|{13:^12}'
    tl1 = ['PDB code','Reported in PDB','Calc from NCS',
           'ASU initial','ASU final','Res.','NCS',
           'Solvent','Data','Year',
           'Use','Geo.','Trans.','time (sec)']
    title1 = t1.format(*tl1)
    # second title row
    t2 = s.replace('.4f','')
    t2 = t2.replace('.2f','')
    t2 = t2.replace(' ','#',1)
    tl2 = [''] + ['r-work','r-free']*4 + ['','copies','fraction']
    tl2 += ['completeness','','NCS','Rest','refine','']
    title2 = t2.format(*tl2)
    summery = title1 + '\n' + title2 + sep + data_str
    # Transformation information
    trans = [list(x.elems) for x in self.init_rotation_matrices]
    rot = [list(x.elems) for x in self.init_translation_vectors]
    transform_before_refinement = [x+y for (x,y) in zip(trans,rot)]
    trans = [x.elems for x in self.rotation_matrices]
    rot = [x.elems for x in self.translation_vectors]
    transform_after_refinement = [x+y for (x,y) in zip(trans,rot)]
    if transform_before_refinement:
      summery_2 = []
      s = '  {0:<37}|{1:^8.4f}|{2:^8.4f}|{3:^8.4f}|{4:^8.4f}|{5:^8.4f}'
      s += '|{6:^8.4f}|{7:^8.4f}|{8:^8.4f}|{9:^8.4f}|{10:^8.4f}|{11:^8.4f}'
      s += '|{12:^8.4f}'
      title1 = ['Rotation,Translation difference','Rxx','Rxy','Rxz','Ryx','Ryy','Ryz','Rzx','Rzy']
      title1.extend(['Rzz','Tx','Ty','Tz'])
      for i,(before, after) in enumerate(zip(transform_before_refinement,
                                      transform_after_refinement)):
        # print actual data
        # t1 = ['Rotation,Translation before'] + list(before)
        # t2 = ['Rotation,Translation after'] + list(after)
        # d1 = s.format(*t1)
        # d2 = s.format(*t2)
        # summery  += d1 + '\n' + d2 + '\n'
        # print difference
        trans_diff = [round(x-y,4) for (x,y) in zip(list(before),list(after))]
        t1 = ['Transform number: {}'.format(i+1)] + trans_diff
        d1 = s.format(*t1)
        if trans_diff.count(0.0) != len(trans_diff):
          summery_2.append(d1)
      if summery_2:
        s = s.replace('.4f','')
        s = s.replace(' ','#',1)
        summery += sep + s.format(*title1) + sep + '\n'.join(summery_2)
    return summery  + '\n'

  def time_from_start(self):
    """ () -> int
    returns the time (in sec) from initiation of refinement object"""
    # consider using self.refinement_time_start instead of self.start_time
    return int(time.time() - self.start_time)

def run(args):
  """
  run tests
  """
  # process args
  help_str = '''\
  To run the test use:
  >>>python test_ncs_refinement.py xxxx
  where xxxx is a 4 characters pdb code

  To run without printing time and refinement steps use:
  >>>python test_ncs_refinement.py xxxx -quiet

  in both cases you will get a table with a summery for the test
  '''
  print_during_refinement = False
  pdb_code = ''
  # make sure will work if calling function with a pdb code string
  if type(args)==str: args = [args]

  if len(args)==0 or ('-h' in args) or ('-help' in args):
    print help_str
    sys.exit()
  elif len(args)>2:
    print help
    raise IOError('To many variables !')
  else:
    for rec in args:
      if len(rec)==4:
        pdb_code = rec
      elif 'quiet' in rec.lower():
        print_during_refinement = False
  assert len(pdb_code) == 4

  # list test to perform
  run_test(pdb_code=pdb_code,
           print_during_refinement=print_during_refinement,
           use_strict_ncs=False,
           use_geometry_restraints = True,
           refine_transformations = False)
  run_test(pdb_code=pdb_code,
           print_during_refinement=print_during_refinement,
           use_strict_ncs=True,
           use_geometry_restraints = True,
           refine_transformations = False)
  run_test(pdb_code=pdb_code,
           print_during_refinement=print_during_refinement,
           use_strict_ncs=True,
           use_geometry_restraints = True,
           refine_transformations = True)


def run_test(pdb_code,
             print_during_refinement,
             use_strict_ncs,
             use_geometry_restraints,
             refine_transformations):
  # Run tests without NCS
  test_obj = ncs_refine_test(
    use_geometry_restraints=use_geometry_restraints,
    real_data=True,
    print_during_refinement=print_during_refinement)
  test_obj.set_working_path()
  test_obj.process_pdb_and_cif_files(pdb_code)
  # Run refinement
  test_obj.refinement_loop(
    n_macro_cycle=5,
    r_work_target=0.1,
    sites=True,
    u_iso=True,
    transformations=refine_transformations,
    use_strict_ncs=use_strict_ncs,
    print_during_refinement=print_during_refinement)
  print test_obj
  del test_obj

if __name__=='__main__':
  # run(sys.argv[1:])

  run('1a37')

  # run('2wws')
  # run('5msf')
  # run('2vq0')

  # Analyze code
  # import cProfile
  # import pstats
  # cProfile.run("run('xxxx')",filename='cProfile_log')
  # p = pstats.Stats('cProfile_log')
  # p.sort_stats('time').print_stats(15)
