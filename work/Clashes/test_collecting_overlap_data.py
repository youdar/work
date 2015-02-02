from __future__ import division
from collecting_overlap_data import Collecting_overlaps_data
import unittest
import sys
import os


csv_data = '''\
PDB ID,Macro mol. overlaps,Sym. overlaps,All overlaps,Macro mol. overlaps per 1000 atoms,Sym. overlaps per 1000 atoms,All overlaps per 1000 atoms,year.,experiment type
# some remarks
1a18,29,5,35,13.9,2.4,16.5,1997,X-RAY DIFFRACTION
xxxx,9,15,85,17.22,0,0,2001,NMR
yyyy,-1,-1,-1,-1,2.4,16.5,2010,X-RAY DIFFRACTION
'''

class Collect_data(unittest.TestCase):

  def setUp(self):
    """ set up paths """
    self.curret_path = os.getcwd()
    self.file_to_delete = []
    if sys.platform.startswith('win'):
      self.test_folder = r'C:\Phenix\Dev\Work\work\Clashes\junk'
    else:
      self.test_folder = '/net/cci-filer2/raid1/home/youval/Work/work/Clashes/junk'
    os.chdir(self.test_folder)


  def test_read_csv_file(self):
    print sys._getframe().f_code.co_name
    fn = 'test_file_csv_1.txt'
    self.file_to_delete.append(fn)
    open(fn,'w').write(csv_data)
    obj = Collecting_overlaps_data()
    obj.data_file_name = os.path.join(self.test_folder,fn)
    self.assertEqual(len(obj.read_csv_data()),3)

  def test_writing_data_files(self):
    print sys._getframe().f_code.co_name
    # create log files
    data = ['1a18,29,5,35,13.9,2.4,16.5,1997,X-RAY DIFFRACTION',
            'xxxx,9,15,85,17.22,0,0,2001,NMR',
            'yyyy,-1,-1,-1,-1,2.4,16.5,2010,X-RAY DIFFRACTION']
    for d in data:
      fn = 'log_' + d[:4] + '.pdb'
      open(fn,'w').write(d)
      self.file_to_delete.append(fn)
    #
    obj = Collecting_overlaps_data()
    self.assertFalse(os.path.isfile(obj.clean_dict_file_name))
    self.assertFalse(os.path.isfile(obj.clean_data_file_name))
    obj.queue_data_path = self.test_folder
    obj.get_test_data()
    self.file_to_delete.append('test_data.txt')
    self.file_to_delete.append('test_clean_dict')
    self.file_to_delete.append('test_clean_data')
    self.assertTrue(os.path.isfile(obj.clean_dict_file_name))
    self.assertTrue(os.path.isfile(obj.clean_data_file_name))
    #
    obj2 = Collecting_overlaps_data()
    obj2.queue_data_path = self.test_folder
    obj2.get_test_data()
    self.assertTrue(len(obj2.data),3)
    self.assertTrue(len(obj2.clean_data),2)


  def tearDown(self):
    """ delete files created in during testing"""
    if self.file_to_delete:
      for fn in self.file_to_delete:
        if os.path.isfile(fn): os.remove(fn)
    os.chdir(self.curret_path)

def run_selected_tests():
  """  Run selected tests

  1) List in "tests" the names of the particular test you want to run
  2) Comment out unittest.main()
  3) Un-comment unittest.TextTestRunner().run(run_selected_tests())
  """
  tests = ['test_writing_data_files']
  suite = unittest.TestSuite(map(Collect_data, tests))
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
