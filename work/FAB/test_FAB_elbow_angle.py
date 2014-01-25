from __future__ import division
from FAB.FAB_elbow_angle import FAB_elbow_angle
from libtbx.utils import null_out
from iotbx.pdb import fetch
import unittest
import cProfile
import shutil
import tempfile
import sys,os
import time

'''
Test Fragment antigen-binding (Fab) elbow angle calcuation

@author Youval Dar (LBL 2014)
'''


class TestFabElbowAngle(unittest.TestCase):

  def setUp(self):
    '''
    create a temporary folder with files for testing.
    '''
    self.currnet_dir = os.getcwd()
    self.tempdir = tempfile.mkdtemp('tempdir')
    os.chdir(self.tempdir)

  def test_1bbd(self):
    '''
    Compare to published value
    http://www.ncbi.nlm.nih.gov/pubmed/16497332
    http://proteinmodel.org/AS2TS/RBOW/calculation.htm
    '''
    fn = '1bbd'
    pdb_fn = fetch.get_pdb (fn,data_type='pdb',mirror='rcsb',log=sys.stdout)

    fab = FAB_elbow_angle(
      pdb_file_name=fn,
      chain_ID_light='L',
      chain_ID_heavy='H',
      limit_light=107,
      limit_heavy=113)

    calculated = fab.FAB_elbow_angle
    expected = 127
    msg = 'FAB angle for {0} is {1:3.0f} instead of {2}'.format(fn,calculated,expected)
    self.assertAlmostEqual(calculated,expected,1, msg)


  def tearDown(self):
    '''remove temp files and folder'''
    #os.chdir(self.currnet_dir)
    #shutil.rmtree(self.tempdir)


if __name__ == "__main__":
  unittest.main(verbosity=2)  # provides a command-line interface to the test script
  #unittest.main()
