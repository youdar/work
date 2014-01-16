from __future__ import division
from iotbx import pdb
import unittest
import cProfile
import shutil
import zipfile
import tempfile
import sys,os

'''
Test strict Non-crystallographic symmetry (NCS)
application in refinement

@author Youval Dar (LBL 2014)
'''


class TestStrictNCS(unittest.TestCase):

  #@classmethod
  #def setUpClass(cls):
    #print 'setUpClass'
    #print '*'*60

  def setUp(self):
    '''
    create a temporary folder with files for testing.
    Using ncs0_pdb, which contains a single NCS, CRYST1 and MTRIX records,
    produce:
    1) asu0.pdb, with complete Asymmetric Unit (ASU)
    2) ncs1.pdb  pertubed version of ncs0_pdb, with a single NCS and MTRIX info
    '''
    self.currnet_dir = os.getcwd()
    self.testdir = tempfile.mkdtemp('testdir')
    os.chdir(self.testdir)
    self.asu0_filename = 'asu0.pdb'
    self.ncs1_filename = 'ncs1.pdb'
    # Create and write a file ncs0.pdb

    # Create and write asu0.pdb
    f = open(self.asu0_filename,'w')
    f.write(ncs0_pdb)
    f.close()
    # Create and write ncs0.pdb
    f = open(self.asu0_filename,'w')
    f.write(ncs0_pdb)
    f.close()

  def test_process_integrity(self):
    ''' Test that when comparing asu0 to itself we get R-work is zero'''
    pdb_inp = pdb.input(file_name=self.asu0_filename)
    pass

  def test_pertubed_ncs(self):
    '''Test that the pertubed NCS (ncs1.pd) is different than the original one (ncs0_pdb)
    by checking that R-work is not zero
    Compare:
    ncs0_pdb to ncs1.pd
    asu0.pdb to asu1.pdb'''
    pass

  def test_refinement(self):
    '''Test that refining asu1.pdb converge to asu0.pdb
    Use asu0.pdb to create x-ray structure (xrs)
    and from it produce F_obs_test (observed structure factors)'''
    pass


  #def test_2(self):
    #'''Test that Error is raised when trying to archive the current directory'''
    ##self.assertRaises(ValueError,archive_dir.archive_dir_files,self.testdir)
    #pass


  def tearDown(self):
    '''remove temp files and folder'''
    os.chdir(self.currnet_dir)
    shutil.rmtree(self.testdir)


# Raw data used to build test cases
ncs0_pdb="""\n
CRYST1   20.000   20.000   20.000  90.00  90.00  90.00 P 1
ATOM      1  N   LYS     1       5.000   5.000   5.000  1.00 20.00           N
ATOM      1  N   LYS     2       6.000   5.000   5.000  1.00 20.00           N
ATOM      1  N   LYS     4       5.000   5.500   5.500  1.00 20.00           N
TER
END
"""

if __name__ == "__main__":
  unittest.main(verbosity=2)  # provides a command-line interface to the test script
  #unittest.main()