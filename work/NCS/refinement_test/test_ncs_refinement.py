from __future__ import division
from iotbx.pdb.multimer_reconstruction import multimer
import mmtbx.refinement.minimization_ncs_constraints
import mmtbx.monomer_library.pdb_interpretation
from libtbx.test_utils import approx_equal
from scitbx.array_family import flex
import mmtbx.monomer_library.server
from libtbx import adopt_init_args
from mmtbx import monomer_library
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
    adopt_init_args(self, locals()) # creates self.xxx for all arguments

  def read_pdb_file(self, file_name):
    """
    Read pdb file
    Check if it contains a single NCS copy
    Do geometry minimization
    """
    os.path.isfile(file_name)

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

def set_working_path():
  """
  Change working directory to avoid littering
  """
  username = getpass.getuser()
  if username.lower() == 'youval':
    osType = sys.platform
    if osType.startswith('win'):
      tempdir = (r'C:\Phenix\Dev\Work\work\NCS\junk\pdb_test')
    else:
      tempdir = ('/net/cci/youval/Work/work/NCS/junk/pdb_test')
    os.chdir(tempdir)


def run (args, out=sys.stdout):
  """
  run test
  """

if __name__=='__main__':
  set_working_path()
  run(sys.argv[1:])