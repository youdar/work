from __future__ import division
import nonbonded_pdb_data_collect as nbo
from libtbx.utils import null_out
from iotbx.pdb import fetch
import unittest
import sys
import os


class Collect_nonbonded_overlap_data(unittest.TestCase):

  def setUp(self):
    """ set up paths """
    self.curret_path = os.getcwd()
    self.file_to_delete = []
    if sys.platform.startswith('win'):
      dr = r'C:\Phenix\Dev\Work\work\Clashes\junk'
    else:
      dr = '/net/cci-filer2/raid1/home/youval/Work/work/Clashes/junk'
    os.chdir(dr)


  def test_get_pdb_id(self):
    """ check that we get clean pdb id code"""
    print sys._getframe().f_code.co_name
    fn1 = r'C:\Phenix\Dev\Work\work\Clashes\junk\2H77.updated.pdb'
    fn2 = 'C:\\Phenix\\Dev\\Work\\work\\Clashes\\junk\\2H77.updated.pdb'
    fn3 = '/net/cci-filer2/raid1/home/youval/Work/work/Clashes/2H77.updated.pdb'
    fn4 = 'C:\\Phenix\\Dev\\Work\\work\\Clashes\\junk\\2H77.pdb'
    fn5 = 'C:\\Phenix\\Dev\\Work\\work\\Clashes\\junk\\2H77'
    fn6 = '2H77.pdb'
    fn7 = '2H77'
    fn8 = 'pdb2H77.ent.gz'
    test_str = [fn1,fn2,fn3,fn4,fn5,fn6,fn7,fn8]
    for fn in test_str:
      file_name = nbo.get_pdb_id(fn)
      self.assertEqual(file_name,'2h77')

  def test_run_good_structure(self):
    """   """
    print sys._getframe().f_code.co_name
    # make sure that the function imports the file
    if os.path.isfile('1a18.pdb'):
      os.remove('1a18.pdb')
    nob_out =  nbo.run('1a18.pdb')
    self.assertEqual(nob_out,'1a18,29,5,35,13.9,2.4,16.5,1997,X-RAY DIFFRACTION')
    nob_out =  nbo.run('1a18')
    self.assertEqual(nob_out,'1a18,29,5,35,13.9,2.4,16.5,1997,X-RAY DIFFRACTION')
    nob_out =  nbo.run('this_is_junk')
    self.assertEqual(nob_out[:21],'this_is_junk,-5,-5,-5')
    self.assertEqual(len(nob_out.split(',')),9)
    # test verbose output
    nob_out = nbo.run('1a18.pdb',verbose=True)
    data = nob_out.split('\n')
    self.assertEqual(len(data),5)
    data = data[-2].split('|')
    data = [x.strip() for x in data]
    data = ','.join(data)
    self.assertEqual(data,'1a18,29,5,35,1997,X-RAY DIFFRACTION')
    self.file_to_delete.append('1a18.pdb')

  def test_catch_unknown_pairs(self):
    print sys._getframe().f_code.co_name
    pdb_id = '3a3w'
    nob_out =  nbo.run(pdb_id)
    expected = [pdb_id] + ['-2']*3
    expected = ','.join(expected)
    self.assertEqual(nob_out[:len(expected)],expected)

  def test_catch_multiple_models(self):
    print sys._getframe().f_code.co_name
    pdb_id = '203d'
    nob_out =  nbo.run(pdb_id)
    expected = [pdb_id] + ['-3']*3
    expected = ','.join(expected)
    self.assertEqual(nob_out[:len(expected)],expected)

  def test_catch_bad_cryst1(self):
    print sys._getframe().f_code.co_name
    pdb_id = '2bvb'
    nob_out =  nbo.run(pdb_id)
    expected = [pdb_id] + ['-4']*3
    expected = ','.join(expected)
    self.assertEqual(nob_out[:len(expected)],expected)

  def tearDown(self):
    """ delete files created in during testing"""
    if self.file_to_delete:
      for fn in self.file_to_delete:
        if os.path.isfile(fn): os.remove(fn)

def run_selected_tests():
  """  Run selected tests

  1) List in "tests" the names of the particular test you want to run
  2) Comment out unittest.main()
  3) Un-comment unittest.TextTestRunner().run(run_selected_tests())
  """
  tests = ['test_catch_bad_cryst1']
  suite = unittest.TestSuite(map(Collect_nonbonded_overlap_data, tests))
  return suite


if __name__ == '__main__':
  try:
    pdb_dir = os.environ["PDB_MIRROR_PDB"]
  except KeyError:
    pdb_dir = ''
  # use for individual tests
  # unittest.TextTestRunner().run(run_selected_tests())

  # Use to run all tests
  unittest.main()
