from __future__ import division
from iotbx.pdb.multimer_reconstruction import multimer
from mmtbx.utils import flex
from cctbx import crystal
from numpy import random
from iotbx import pdb
import shutil
import tempfile
import sys,os

'''
Create test file for testing strict Non-crystallographic symmetry (NCS)
application in refinement

@author Youval Dar (LBL 2014)
'''


def run(ncs0_pdb):
  '''
  create a temporary folder with files for testing.
  Using ncs0_pdb, which contains a single NCS, CRYST1 and MTRIX records,
  produce:
  1) ncs0.pdb, with CRYST1 and SCALE records
  2) asu0.pdb, with complete Asymmetric Unit (ASU)
  3) ncs1.pdb  pertubed version of ncs0_pdb, with a single NCS and MTRIX info

  Then print all those files on screen, so that they can be copy and paste
  to the test_strict_NCS.py file
  '''
  currnet_dir = os.getcwd()
  tempdir = tempfile.mkdtemp('tempdir')
  os.chdir(tempdir)
  #
  testdir = r'C:\Phenix\Dev\Work\work\NCS\test_files'
  os.chdir(testdir)
  #
  print 'working at {}\n'.format(os.getcwd())
  #
  asu0_filename = 'asu0.pdb'
  ncs1_filename = 'ncs1.pdb'
  ncs0_filename = 'ncs0.pdb'
  asu1_filename = 'asu1.pdb'
  #
  f = open(ncs0_filename,'w').write(ncs0_pdb)
  f = open('ncs0_origin.pdb','w').write(ncs0_pdb)
  # Create the ASU coordinates using MTRIX records
  crystal_symmetry = create_asu(ncs_filename=ncs0_filename, asu_filename=asu0_filename)
  # Add CRYST1 records to ncs0
  pdb_inp = pdb.input(file_name = ncs0_filename)

  #xrs = pdb_inp.xray_structure_simple()
  #crystal_symmetry = xrs.crystal_symmetry()

  pdb_inp.write_pdb_file(file_name = ncs0_filename, crystal_symmetry = crystal_symmetry)
  # When using pdb_inp.write_pdb_file the MTRIX record are omitted. Add them back
  add_MTRIX_to_pdb(ncs0_filename, 'ncs0_origin.pdb')
  print '*'*100
  print 'The generating NCS'
  print_file(ncs0_filename)
  print '*'*100
  print 'Complete Asymmetric Unit (ASU) - use as target'
  print_file(asu0_filename)
  print 'ncs0 need to have the same CRYST1 and SCALE records and asu0'
  print '*'*100
  #
  # Shake ncs0 to create ncs1
  pdb_inp = pdb.input(file_name = ncs0_filename)
  xyz = pdb_inp.atoms().extract_xyz()
  sig = 0.5
  # apply random shift on all atoms
  xyz = random.normal(xyz,sig)
  tmp = [tuple(x) for x in xyz]
  pdb_inp.atoms().set_xyz(flex.vec3_double(tmp))
  pdb_inp.write_pdb_file(file_name = ncs1_filename, crystal_symmetry = crystal_symmetry)
  print 'Shaked coordinates'
  print '-'*100
  print 'Use the coordinates below for the modified NCS'
  # add MTRIX records to file
  add_MTRIX_to_pdb(ncs1_filename, ncs0_filename)
  print_file(ncs1_filename)

  # Create the ASU from the modified coordinates
  create_asu(ncs_filename=ncs1_filename, asu_filename=asu1_filename, crystal_symmetry=crystal_symmetry)
  print '*'*100

  print 'Done.'
  # Cleanup
  os.chdir(currnet_dir)
  shutil.rmtree(tempdir)

def add_MTRIX_to_pdb(pdb_fn, pdb_MTRIX_fn):
  '''
  Add MTRIX records from the file pdb_MTRIX_fn
  to the file pdb_fn and write the modified pdb file, with
  the MTRIX records, in the current directory
  '''
  ncs1_pdb = open(pdb_fn,'r').read().splitlines()
  ncs0_pdb = open(pdb_MTRIX_fn,'r').read().splitlines()
  i = 0
  while ncs0_pdb:
    x = ncs0_pdb.pop(0)
    if x.startswith('MTRIX'):
      ncs1_pdb.insert(i+4, x)
      i += 1
  ncs1_pdb = '\n'.join(ncs1_pdb)
  f = open(pdb_fn,'w').write(ncs1_pdb)

def create_asu(ncs_filename,asu_filename,crystal_symmetry=None):
  ''' (str,str) -> crystal_symmetry object
  Create ASU from NCS and save the new pdb file
  with CRYST1 and SCALE records in local directory

  Argument:
  ---------
  ncs_filename : NCS to be used to create the ASU
  asu_filename : ASU file name
  crystal_symmetry : Allow forcing other crystal_symmetry

  Returns:
  --------
  crystal_symmetry
  '''
  m = multimer(ncs_filename,'cau',error_handle=True,eps=1e-2)
  if m.number_of_transforms == 0:
    print 'Number of transforms is zero'
  m.write(asu_filename)
  pdb_inp = pdb.input(file_name = asu_filename)
  xrs = pdb_inp.xray_structure_simple()
  if not crystal_symmetry:
    crystal_symmetry = xrs.crystal_symmetry()
  pdb_inp.write_pdb_file(file_name = asu_filename, crystal_symmetry = crystal_symmetry)
  return crystal_symmetry


def print_file(fn):
  print '='*100
  print open(fn,'r').read()
  print '='*100

if __name__ == "__main__":

  # Raw data used to build test cases
  # The coordinates below where first optimized using phenix.geometry_minimization
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
ATOM      1  N   THR 1   1       9.670  10.289  11.135  1.00 26.11           N
ATOM      2  CA  THR 1   1       9.559   8.931  10.615  1.00 27.16           C
ATOM      3  C   THR 1   1       9.634   7.903  11.739  1.00 20.29           C
ATOM      4  O   THR 1   1      10.449   8.027  12.653  1.00 35.00           O
ATOM      5  CB  THR 1   1      10.660   8.630   9.582  1.00 34.84           C
ATOM      6  OG1 THR 1   1      10.560   9.552   8.490  1.00 67.35           O
ATOM      7  CG2 THR 1   1      10.523   7.209   9.055  1.00 43.17           C
TER
"""

  run(ncs0_pdb)