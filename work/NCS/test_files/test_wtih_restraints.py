from __future__ import division
import mmtbx.monomer_library.pdb_interpretation
from mmtbx import monomer_library
import mmtbx.monomer_library.server
import getpass
import sys
import os

def run(file_name):
  pdb_processed_file = monomer_library.pdb_interpretation.run(args=[file_name],
    assume_hydrogens_all_missing=False,
    hard_minimum_nonbonded_distance=0.0,
    nonbonded_distance_threshold=None,
    substitute_non_crystallographic_unit_cell_if_necessary=True)

  grm = pdb_processed_file.geometry_restraints_manager()
  print 'done'

def set_test_folder():
    """
    Change working directory to avoid littering of
    phenix_sources\phenix_regression\development\ncs_constraints.py
    """
    username = getpass.getuser()
    if username.lower() == 'youval':
      osType = sys.platform
      if osType.startswith('win'):
        tempdir = (r'C:\Phenix\Dev\Work\work\NCS\junk')
      else:
        tempdir = ('/net/cci/youval/Work/work/NCS/junk')
      os.chdir(tempdir)

if __name__=='__main__':
  set_test_folder()
  run('full_asu.pdb')