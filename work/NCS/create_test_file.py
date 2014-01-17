from __future__ import division
from iotbx.pdb.multimer_reconstruction import multimer
from mmtbx import utils
from iotbx import reflection_file_reader
from libtbx import easy_run
from iotbx import pdb
import unittest
import cProfile
import shutil
import zipfile
import tempfile
from cctbx import crystal
import sys,os

'''
Create test file for testing strict Non-crystallographic symmetry (NCS)
application in refinement

@author Youval Dar (LBL 2014)
'''

def run():
  '''
  create a temporary folder with files for testing.
  Using ncs0_pdb, which contains a single NCS, CRYST1 and MTRIX records,
  produce:
  1) ncs0.pdb, with CRYST1 and SCALE records
  2) asu0.pdb, with complete Asymmetric Unit (ASU)
  3) ncs1.pdb  pertubed version of ncs0_pdb, with a single NCS and MTRIX info
  '''
  currnet_dir = os.getcwd()
  tempdir = tempfile.mkdtemp('tempdir')
  os.chdir(tempdir)
  #
  testdir = r'C:\Phenix\Dev\Work\work\NCS\test_files'
  os.chdir(testdir)
  #
  print 'working at {}'.format(os.getcwd())
  #
  asu0_filename = 'asu0.pdb'
  ncs1_filename = 'ncs1.pdb'
  ncs0_filename = 'ncs0.pdb'
  #
  f = open(ncs0_filename,'w')
  f.write(ncs0_pdb)
  f.close()
  # Create the ASU coordinates using MTRIX records
  # Do it before creating the CRYST1 records, when creating them the MTRIX will not be saved
  m = multimer(ncs0_filename,'cau',error_handle=True,eps=1e-2)
  if m.number_of_transforms == 0:
    print 'Number of transforms is zero'
  m.write(asu0_filename)
  pdb_inp = pdb.input(file_name = asu0_filename)
  xrs = pdb_inp.xray_structure_simple()
  crystal_symmetry = xrs.crystal_symmetry()
  pdb_inp.write_pdb_file(file_name = asu0_filename, crystal_symmetry = crystal_symmetry)
  pdb_inp.crystal_symmetry()
  print '*'*100
  print 'Complete Asymmetric Unit (ASU) - use as target'
  get_file_as_str(asu0_filename)
  # Create and write a file ncs0.pdb
  pdb_inp = pdb.input(file_name = ncs0_filename)
  xrs = pdb_inp.xray_structure_simple()
  crystal_symmetry = xrs.crystal_symmetry()
  pdb_inp.write_pdb_file(file_name = ncs0_filename, crystal_symmetry = crystal_symmetry)
  print '*'*100
  print 'For the NCS remember to add the MTRIX portion (see below) back to the pdb'
  get_file_as_str(ncs0_filename)
  print 'The MTRIX records for the NCS'
  print ncs0_pdb
  print 'Done.'

  # Cleanup
  os.chdir(currnet_dir)
  shutil.rmtree(tempdir)

def get_file_as_str(fn):
  print '='*100
  print open(fn,'r').read()
  print '='*100



# Raw data used to build test cases
ncs0_pdb="""\
MTRIX1   1  1.000000  0.000000  0.000000        0.00000    1
MTRIX2   1  0.000000  1.000000  0.000000        0.00000    1
MTRIX3   1  0.000000  0.000000  1.000000        0.00000    1
MTRIX1   2  0.496590 -0.643597  0.582393        0.00000
MTRIX2   2  0.867925  0.376088 -0.324443        0.00000
MTRIX3   2 -0.010221  0.666588  0.745356        0.00000
MTRIX1   3 -0.317946 -0.173437  0.932111        0.00000
MTRIX2   3  0.760735 -0.633422  0.141629        0.00000
MTRIX3   3  0.565855  0.754120  0.333333        0.00000
ATOM      1  N   THR 1   1       9.483  10.256  10.995  1.00 26.11           N
ATOM      2  CA  THR 1   1       9.489   8.820  10.603  1.00 27.16           C
ATOM      3  C   THR 1   1       9.725   7.911  11.796  1.00 20.29           C
ATOM      4  O   THR 1   1      10.586   8.189  12.629  1.00 35.00           O
ATOM      5  CB  THR 1   1      10.632   8.506   9.642  1.00 34.84           C
ATOM      6  OG1 THR 1   1      10.831   9.604   8.746  1.00 67.35           O
ATOM      7  CG2 THR 1   1      10.308   7.254   8.859  1.00 43.17           C
TER
"""

if __name__ == "__main__":
  run()