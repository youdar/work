from __future__ import division
from iotbx.pdb.multimer_reconstruction import multimer
import mmtbx.refinement.minimization_ncs_constraints
import mmtbx.monomer_library.pdb_interpretation
from libtbx.test_utils import approx_equal
from scitbx.array_family import flex
import mmtbx.monomer_library.server
from libtbx import adopt_init_args
from libtbx.utils import null_out
from mmtbx import monomer_library
from iotbx.pdb import fetch
from cctbx import xray
import mmtbx.f_model
import mmtbx.utils
import iotbx.cif
import iotbx.pdb
import getpass
import os
import sys

__author__ = 'Youval Dar'

class ncs_refine_test(object):
  def __init__(self,
               n_macro_cycle = 1,
               sites = True,
               u_iso = False,
               alternate_refinement = True,
               finite_grad_differences_test = False,
               use_geometry_restraints = True,
               r_work_target = 0.01,
               d_min = 0.2):
    """
    Arguments:
    n_macro_cycle: (int) number of refinement cycles
    sites: (bool) refinement using sites
    u_iso: (bool) refinement using b-factors (ADP)
    alternate_refinement: (bool) when True, each refinement cycle will run
                          site refinement and u_iso (ADP) refinement
    finite_grad_differences_test: (bool) Compare calculated gradient to
                                  estimated one
    use_geometry_restraints: (bool) Use geometry restraints in refinement
    r_work_target: (float) stop refinement when r_work <= r_work_target
    d_min: (float)

    """
    # creates self.xxx for all arguments
    adopt_init_args(self, locals())
    # Check for invalid input parameters
    if [sites, u_iso].count(True) != 1:
      raise IOError('Refinement method error')
    if n_macro_cycle > 100:
      raise IOError('To many refinement cycles')
    # check if working at LBL local machine
    osType = sys.platform
    self.pdb_dir = ''
    if osType.startswith('win'):
      self.work_env = 'win'
    else:
      try:
        self.pdb_dir = os.environ["PDB_MIRROR_PDB"]
        self.cif_dir = os.environ["PDB_MIRROR_STRUCTURE_FACTORS"]
        self.work_env = 'lbl'
      except KeyError:
        self.work_env = 'other'

  def process_pdb_and_cif_files(self, args):
    """
    Read pdb file when not working at LBL

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
    else:
      self.full_path_pdb = self.get_full_path(data_type='pdb')
    # Get cif file
    if os.path.isfile(self.cif_file_name):
      self.full_path_cif = os.path.realpath(self.cif_file_name)
    else:
      self.full_path_cif = self.get_full_path(data_type='xray')
    # Collect Fobs and r_free_flags from cif file
    self.get_structure_factors()
    # Process PDB file
    self.process_pdb_file()

  def process_pdb_file(self, build_unit_cell=False):
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
    self.is_single_ncs = (m.number_of_transforms > 0)
    #
    pdb_inp = iotbx.pdb.input(file_name=self.full_path_pdb)
    ph_ncs = pdb_inp.construct_hierarchy()
    xrs_ncs = pdb_inp.xray_structure_simple()
    crystal_symmetry = xrs_ncs.crystal_symmetry()
    # Collect r-work for the single NCS copy
    self.get_pdb_ncs_r_work_values(pdb_inp=pdb_inp,xrs=xrs_ncs)
    xrs_asu = m.assembled_multimer.extract_xray_structure(
      crystal_symmetry = crystal_symmetry)
    # If needed, construct unit cell
    if build_unit_cell:
      xrs_unit_cell = xrs_asu.orthorhombic_unit_cell_around_centered_scatterers(
        buffer_size=2)
      crystal_symmetry = xrs_unit_cell.crystal_symmetry()
      xrs_ncs = pdb_inp.xray_structure_simple(
        crystal_symmetry=crystal_symmetry)
      xrs_asu = m.assembled_multimer.extract_xray_structure(
        crystal_symmetry = crystal_symmetry)
      ph_ncs.adopt_xray_structure(xrs_ncs)
    assert xrs_ncs.crystal_symmetry().\
      is_similar_symmetry(xrs_asu.crystal_symmetry())
    self.ncs_hierarchy = ph_ncs
    self.crystal_symmetry = crystal_symmetry
    self.rotation_matrices = m.rotation_matrices
    self.translation_vectors = m.translation_vectors
    selection_str = 'chain ' + ' and chain '.join(m.ncs_chains_ids)
    self.ncs_selection = m.assembled_multimer.\
      atom_selection_cache().selection(selection_str)
    # Get geometry restraints manager (grm)
    self.grm = None
    pdb_str = m.assembled_multimer.as_pdb_string(
      crystal_symmetry=crystal_symmetry)
    if self.use_geometry_restraints:
      self.grm = get_restraints_manager(pdb_string=pdb_str)
    # Get f_model and initial r_work
    self.fmodel, self.initial_r_work = self.get_f_model(xrs=xrs_asu)

  def get_pdb_ncs_r_work_values(self,pdb_inp, xrs):
    """
    Get the r-work value reported in PDB
    Calculate r-work from a single NCS copy (from f-model)

    Arguments:
    pdb_inp: PDB input object of the single NCS
    xrs: x-ray structure object for a single NCE
    """
    pio = pdb_inp.get_r_rfree_sigma()
    self.r_work_reported_pdb_ncs = pio.r_work
    _,self.r_work_calc_pdb_ncs = self.get_f_model(xrs)

  def get_structure_factors(self):
    """
    Get f_obs and r_free_flags
    From cif file if available
    """
    self.f_obs = self.r_free_flags = None
    if self.full_path_cif:
      miller_arrays = iotbx.cif.reader(file_path=self.full_path_cif).as_miller_arrays()
      for ma in miller_arrays:
        # Consider using Bijvoet mates
        ma = ma.average_bijvoet_mates()
        if ma.is_xray_amplitude_array():
          self.f_obs = abs(ma)
          break
        elif ma.is_xray_intensity_array():
          # convert i_obs to f_obs
          self.f_obs = abs(ma.french_wilson(log=null_out()))
    else:
      raise RuntimeError("No cif file.")

    if self.f_obs:
      self.r_free_flags = self.f_obs.generate_r_free_flags()
      # self.f_obs.show_summary()
      # generate mtz file - Un-comment if needed
      # mtz_dataset = self.f_obs.as_mtz_dataset(column_root_label="F-obs")
      # mtz_dataset.add_miller_array(
      #   miller_array=self.r_free_flags,
      #   column_root_label="R-free-flags")
      # mtz_object = mtz_dataset.mtz_object()
      # mtz_object.write(file_name = 'data_{}.mtz'.format(self.pdb_code))
      # del mtz_object,mtz_dataset
    else:
      raise RuntimeError("Missing amplitude array.")

  def get_f_model(self, xrs):
    """
    Argument:
    xrs: x-ray structure factor object

    Return: fmodel,r_work
    fmodel: f-model object
    r_work: (float) r-work derived from the f-model
    """
    params = mmtbx.f_model.sf_and_grads_accuracy_master_params.extract()
    params.algorithm = "direct"
    fmodel = mmtbx.f_model.manager(
      f_obs                        = self.f_obs,
      r_free_flags                 = self.r_free_flags,
      xray_structure               = xrs,
      sf_and_grads_accuracy_params = params,
      target_name                  = "ls_wunit_k1")
    # Update model for Bulk-solvent and overall scaling
    fmodel.update_all_scales(update_f_part1_for=None)
    r_work = fmodel.r_work()
    return fmodel, r_work

  def refine_using_complete_asu(self):
    """
    Refine the complete ASU using strict-NCS refinement
    """
    if self.alternate_refinement:
      self.n_macro_cycle *= 2
    for macro_cycle in xrange(self.n_macro_cycle):
      if self.alternate_refinement and macro_cycle%2==1:
        # alternate refinement type
        self.sites = not self.sites
        self.u_iso = not self.u_iso
      data_weight = None
      if(self.use_geometry_restraints):
        data_weight = self.get_weight()
      minimized = mmtbx.refinement.minimization_ncs_constraints.lbfgs(
        fmodel                       = self.fmodel,
        rotation_matrices            = self.rotation_matrices,
        translation_vectors          = self.translation_vectors,
        ncs_atom_selection           = self.ncs_selection,
        finite_grad_differences_test = self.finite_grad_differences_test,
        geometry_restraints_manager  = self.grm,
        data_weight                  = self.get_weight(),
        refine_sites                 = self.sites,
        refine_u_iso                 = self.u_iso)
      refine_type = 'adp'*self.u_iso + 'sites'*self.sites
      print "  macro_cycle %3d (%s)   r_factor: %6.4f   " \
            "finite_grad_difference_val: %.4f"%(
        macro_cycle, refine_type,
        self.fmodel.r_work(),
        minimized.finite_grad_difference_val)
      if self.fmodel.r_work()<=self.r_work_target: break

    self.fmodel = self.get_f_model()
  def get_weight(self):
    fmdc = self.fmodel.deep_copy()
    fmdc.xray_structure.shake_sites_in_place(mean_distance=0.3)
    fmdc.update_xray_structure(xray_structure = fmdc.xray_structure,
      update_f_calc=True)
    fmdc.xray_structure.scatterers().flags_set_grads(state=False)
    xray.set_scatterer_grad_flags(
      scatterers = fmdc.xray_structure.scatterers(),
      site       = True)
    # fmodel gradients
    gxc = flex.vec3_double(fmdc.one_time_gradients_wrt_atomic_parameters(
      site = True).packed())
    # manager restraints, energy sites gradients
    gc = self.grm.energies_sites(
      sites_cart        = fmdc.xray_structure.sites_cart(),
      compute_gradients = True).gradients
    gc_norm  = gc.norm()
    gxc_norm = gxc.norm()
    weight = 1.
    if(gxc_norm != 0.0):
      weight = gc_norm / gxc_norm
    return  weight

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
    full_path = ''
    if self.work_env == 'lbl':
      # Use file directly from LBL MIRRORs
      if file_type=='pdb':
        mirror_dir = self.pdb_dir
        s=-12; e=-8  # start and end of 4 letter pdb code in record
      elif file_type=='cif':
        mirror_dir = self.cif_dir
        s=-14; e=-10  # start and end of 4 letter pdb code in record
      else: raise IOError('Wrong file format')
      pdb_files = open(os.path.join(mirror_dir, "INDEX"), "r").readlines()
      # assuming that the 4 letter PDB code is in [s:e] of the records
      files = [x for x in pdb_files if x[s:e]==self.pdb_code]
      if len(files)!=1:
        outstr = 'Found {0} number of {1} in PDB_MIRROR_PDB index'
        print outstr.format(len(files),self.pdb_code)
      else:
        full_path = os.path.join(self.tempdir,files[0]).strip()
    else:
      # Fetch file to local current folder
      fetched_file = fetch.get_pdb (
        id=self.pdb_code, data_type=data_type,
        mirror='rcsb',quiet=True,log=sys.stdout)
      full_path = os.path.realpath(fetched_file)
    return full_path

  def set_working_path(self):
    """
    Change working directory to avoid littering
    """
    self.current_folder = os.getcwd()
    username = getpass.getuser()
    if username.lower() == 'youval':
      osType = sys.platform
      if osType.startswith('win'):
        self.tempdir = (r'C:\Phenix\Dev\Work\work\NCS\junk\pdb_test')
      else:
        self.tempdir = ('/net/cci/youval/Work/work/NCS/junk/pdb_test')
      os.chdir(self.tempdir)

def get_restraints_manager(pdb_file_name=None,pdb_string=None):
  assert [pdb_file_name,pdb_string].count(None)==1
  mon_lib_srv = monomer_library.server.server()
  ener_lib = monomer_library.server.ener_lib()
  if pdb_string: pdb_lines = pdb_string.splitlines()
  else: pdb_lines = None
  processed_pdb_file = monomer_library.pdb_interpretation.process(
    mon_lib_srv    = mon_lib_srv,
    ener_lib       = ener_lib,
    file_name      = pdb_file_name,
    raw_records    = pdb_lines,
    force_symmetry = True)
  geometry = processed_pdb_file.geometry_restraints_manager(
    show_energies = False, plain_pairs_radius = 5.0)
  return mmtbx.restraints.manager(
    geometry = geometry, normalization = False)

def run (args, out=sys.stdout):
  """
  run tests
  """
  test_obj = ncs_refine_test()
  test_obj.set_working_path()
  test_obj.process_pdb_and_cif_files(args)
  outstr = 'r-factors: reported in pdb: {0:<10.5f}from ncs fmodel: {1:<10' \
           '.5f}initial value for ASU fmodel: {2:<10.5f}'
  outstr.format(test_obj.r_work_reported_pdb_ncs, test_obj
                .r_work_calc_pdb_ncs,test_obj.initial_r_work)
  print outstr
  # Run refinement
  test_obj.refine_using_complete_asu()
  print 'Done'

if __name__=='__main__':
  """
  List of files to try

  3p0s, 2wws, 2e0z
  1b35, 1za7, 2buk

  """
  run('3p0s')
  # run(sys.argv[1:])

  # Analyze code
  # import cProfile
  # import pstats
  # cProfile.run("run('3p0s')",filename='cProfile_log')
  # p = pstats.Stats('cProfile_log')
  # p.sort_stats('time').print_stats(15)
