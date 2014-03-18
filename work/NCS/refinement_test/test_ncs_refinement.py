from __future__ import division
from iotbx.pdb.multimer_reconstruction import multimer
import mmtbx.refinement.minimization_ncs_constraints
import mmtbx.monomer_library.pdb_interpretation
from libtbx.test_utils import approx_equal
from scitbx.array_family import flex
import mmtbx.monomer_library.server
from libtbx import adopt_init_args
from mmtbx import monomer_library
from iotbx.pdb import fetch
from cctbx import xray
import mmtbx.f_model
import mmtbx.utils
import iotbx.pdb
import getpass
import os
import sys

__author__ = 'Youval Dar'

class ncs_refine_test(object):
  def __init__(self,
               n_macro_cycle = 10,
               sites = True,
               u_iso = False,
               finite_grad_differences_test = False,
               use_geometry_restraints = True,
               d_min = 0.2):
    """
    Arguments:
    n_macro_cycle: (int) number of refinement cycles
    sites: (bool) refinement using sites
    u_iso: (bool) refinement using b-factors (ADP)
    finite_grad_differences_test: (bool) Compare calculated gradient to
                                  estimated one
    use_geometry_restraints: (bool) Use geometry restraints in refinement
    d_min: (float)

    """
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
        self.work_env = 'lbl'
      except KeyError:
        self.work_env = 'other'




  def read_pdb_file(self, file_name):
    """
    Read pdb file when not working at LBL

    Argument:
    file_name: (str) 4 letters pdb code
    """
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
    # Process PDB file
    self.process_pdb_file()


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


  def process_pdb_file(self, build_unit_cell=False):
    """
    Process file
    - Check if it contains a single NCS copy
    - Build xray structure
    - Get Crystal symmetry

    Arguments:
    build_unit_cell: (bool) leave original crystal symmetry or
                     create a unit cell

    Builds attributes:
    PDB Hierarchy, MTRIX transformation info and NCS atoms selection
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
    xrs_asu = m.assembled_multimer.extract_xray_structure(
      crystal_symmetry = crystal_symmetry)
    # If needed, construct unit cell
    if build_unit_cell:
      xrs_unit_cell = xrs_asu.orthorhombic_unit_cell_around_centered_scatterers(
        buffer_size=2)
      crystal_symmetry = xrs_unit_cell.crystal_symmetry()
      xrs_ncs = pdb_inp.xray_structure_simple(crystal_symmetry=crystal_symmetry)
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

  def process_ncs(self):
    """
    Do geometry minimization to the NCS copy
    Produce ASU
    Create a cell around ASU
    """

  def refine_using_strict_ncs(self):
    """

    """

  def refine_using_complete_asu(self):
    """

    """

  def cleanup(self):
    """
    Cleanup temporary files and folders created during test
    """
    os.chdir(self.current_folder)

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


def run (args, out=sys.stdout):
  """
  run tests
  """
  test_obj = ncs_refine_test()
  test_obj.set_working_path()
  test_obj.read_pdb_file(args)


  print 'Done'

if __name__=='__main__':
  """
  List of files to try

  3p0s, 2wws, 2e0z
  1b35, 1za7, 2buk

  """
  run('3p0s')
  # run(sys.argv[1:])