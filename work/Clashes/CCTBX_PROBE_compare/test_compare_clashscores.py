from __future__ import division
import compare_clashscores_CCTBX_vs_PROBE as cvp
from libtbx.utils import null_out
from iotbx.pdb import fetch
import unittest
import sys
import os

__author__ = 'Youval'


class Compare_clashscores_CCTBX_vs_PROBE(unittest.TestCase):

  def setUp(self):
    """ set up paths """
    self.curret_path = os.getcwd()
    self.file_to_delete = []
    if sys.platform.startswith('win'):
      dr = r'C:\Phenix\Dev\Work\work\Clashes\junk'
    else:
      dr = '/net/cci-filer2/raid1/home/youval/Work/work/Clashes/junk'
    os.chdir(dr)

  def test_get_clashscore_internal(self):
    """ check that we get CCTBX clashscores """
    print sys._getframe().f_code.co_name
    fn = '1G9W'
    # fetch pdb file
    fn = fetch.get_pdb (fn,'pdb',mirror='rcsb',log=null_out())
    self.file_to_delete.append(fn)
    simple,sym,solv,all =  cvp.get_cctbx_clashscores(file_name=fn)
    self.assertEqual(round(simple+sym+solv),round(all))

  def test_get_updated_file_name(self):
    """ check updating file using phenix.reduce """
    print sys._getframe().f_code.co_name
    fn_start = '2H77'
    # fetch pdb file
    fn = fetch.get_pdb (fn_start,'pdb',mirror='rcsb',log=null_out())
    self.file_to_delete.append(fn)
    update_fn = cvp.get_updated_file_name(file_name=fn)
    expected = fn_start + '.updated.pdb'
    self.assertEqual(expected,update_fn)
    self.file_to_delete.append(update_fn)

  def test_get_probe_clashscore(self):
    """ check that we get PROBE clashscores """
    print sys._getframe().f_code.co_name
    fn = '1G9W'
    # fetch pdb file
    fn = fetch.get_pdb (fn,'pdb',mirror='rcsb',log=null_out())
    self.file_to_delete.append(fn)
    clashscore =  cvp.get_probe_clashscore(file_name=fn)
    self.assertEqual(round(clashscore),0)

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
      file_name = cvp.get_pdb_id(fn)
      self.assertEqual(file_name,'2h77')

  def test_get_pdb_file(self):
    """ that the import of file when it does not exist in local folder  """
    print sys._getframe().f_code.co_name
    fn = '3dar.pdb'
    pdb_id = cvp.get_pdb_id(fn)
    have_file = os.path.isfile(fn)
    # Make sure file is not already in local folder
    self.assertFalse(have_file)
    file_name = cvp.get_file_from_rcsb(pdb_id)
    have_file = os.path.isfile(file_name)
    self.assertTrue(have_file)
    self.file_to_delete.append(file_name)

  def test_run_all(self):
    """ Check that we get both CCTBX and PROBE clashscores  """
    print sys._getframe().f_code.co_name
    self.assertFalse(os.path.isfile('1a18.pdb'))
    clashscores =  cvp.run('1a18.pdb')
    self.assertEqual(clashscores,'1a18::13.7::2.4::0.0::16.1::13.0')
    clashscores =  cvp.run('1a18')
    self.assertEqual(clashscores,'1a18::13.7::2.4::0.0::16.1::13.0')
    clashscores =  cvp.run('this_is_junk')
    self.assertEqual(clashscores,'::-1.0::-1.0::-1.0::-1.0::-1.0')
    self.assertEqual(len(clashscores.split('::')),6)
    # test verbose output
    clashscores = cvp.run('1a18.pdb',verbose=True)
    data = clashscores.split('\n')
    self.assertEqual(len(data),6)
    data = data[-2].split()
    data = [x for x in data if x != '|']
    data = '::'.join(data)
    self.assertEqual(data,'1a18::13.7::2.4::0.0::16.1::13.0')
    self.file_to_delete.append('1a18.pdb')

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
  tests = ['test_print']
  suite = unittest.TestSuite(map(Compare_clashscores_CCTBX_vs_PROBE, tests))
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
