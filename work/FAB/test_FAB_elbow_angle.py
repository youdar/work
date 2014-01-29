from __future__ import division
from FAB.FAB_elbow_angle import FAB_elbow_angle
from libtbx.utils import null_out
from scitbx import matrix
from iotbx.pdb import fetch
from math import cos
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

    Tests for 1bbd,7fab,1dba,1plg,1nl0
    are based on:
     Stanfield, et al., JMB 2006
    doi:10.1016/j.jmb.2006.01.023
    http://www.ncbi.nlm.nih.gov/pubmed/16497332
    http://proteinmodel.org/AS2TS/RBOW/calculation.htm
    and
    pymol
    http://www.pymolwiki.org/index.php/Elbow_angle
    '''
    self.currnet_dir = os.getcwd()
    self.tempdir = tempfile.mkdtemp('tempdir')
    os.chdir(self.tempdir)
    # Set delta for testing angles (degrees)
    self.delta = 7

    # Remove this os.chdir when test is working
    os.chdir(r'C:\Phenix\Dev\Work\work\FAB')

  @unittest.skip('Skip test')
  def test_1bbd(self):
    '''Compare to published value'''
    fn = '1bbd'
    pdb_fn = fetch.get_pdb (fn,data_type='pdb',mirror='rcsb',log=sys.stdout)

    fab = FAB_elbow_angle(
      pdb_file_name=fn,
      chain_ID_light='L',
      chain_ID_heavy='H',
      limit_light=114,
      limit_heavy=118)

    calculated = fab.FAB_elbow_angle
    expected = 127
    pymol = 125
    msg = 'FAB angle for {0} is {1:3.0f} instead of {2}'.format(fn,calculated,expected)
    self.assertAlmostEqual(calculated,expected,delta=self.delta,msg=msg)

  @unittest.skip('Skip test')
  def test_7fab(self):
    '''Compare to published value'''
    fn = '7fab'
    pdb_fn = fetch.get_pdb (fn,data_type='pdb',mirror='rcsb',log=sys.stdout)

    fab = FAB_elbow_angle(
      pdb_file_name=fn,
      chain_ID_light='L',
      chain_ID_heavy='H',
      limit_light=104,
      limit_heavy=117)

    calculated = fab.FAB_elbow_angle
    expected = 132
    pymol = 126
    msg = 'FAB angle for {0} is {1:3.0f} instead of {2}'.format(fn,calculated,expected)
    self.assertAlmostEqual(calculated,expected,delta=self.delta,msg=msg)

  @unittest.skip('Skip test')
  def test_1dba(self):
    '''Compare to published value'''
    fn = '1dba'
    pdb_fn = fetch.get_pdb (fn,data_type='pdb',mirror='rcsb',log=sys.stdout)
    fab = FAB_elbow_angle(pdb_file_name=fn)
    calculated = fab.FAB_elbow_angle
    expected = 183
    pymol = 176
    msg = 'FAB angle for {0} is {1:3.0f} instead of {2}'.format(fn,calculated,expected)
    self.assertAlmostEqual(calculated,expected,delta=self.delta,msg=msg)

  @unittest.skip('Skip test')
  def test_1plg(self):
    '''Compare to published value'''
    fn = '1plg'
    pdb_fn = fetch.get_pdb (fn,data_type='pdb',mirror='rcsb',log=sys.stdout)

    fab = FAB_elbow_angle(
      pdb_file_name=fn,
      chain_ID_light='L',
      chain_ID_heavy='H',
      limit_light=112,
      limit_heavy=117)

    calculated = fab.FAB_elbow_angle
    expected = 190
    pymol = 189
    msg = 'FAB angle for {0} is {1:3.0f} instead of {2}'.format(fn,calculated,expected)
    self.assertAlmostEqual(calculated,expected,delta=self.delta,msg=msg)

  #@unittest.skip('Skip test')
  def test_1nl0(self):
    '''Compare to published value'''
    fn = '1nl0'
    pdb_fn = fetch.get_pdb (fn,data_type='pdb',mirror='rcsb',log=sys.stdout)
    fab = FAB_elbow_angle(pdb_file_name=fn)
    calculated = fab.FAB_elbow_angle
    expected = 220
    pymol = 221
    msg = 'FAB angle for {0} is {1:3.0f} instead of {2}'.format(fn,calculated,expected)
    self.assertAlmostEqual(calculated,expected,delta=self.delta,msg=msg)


  def tearDown(self):
    '''remove temp files and folder'''
    #os.chdir(self.currnet_dir)
    #shutil.rmtree(self.tempdir)


if __name__ == "__main__":
  unittest.main(verbosity=2)  # provides a command-line interface to the test script
  #unittest.main()
