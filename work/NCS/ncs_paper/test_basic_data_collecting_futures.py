from __future__ import division
from datetime import datetime
import collect_ncs_files
import unittest
import shutil
import os

__author__ = 'Youval'

# control if temp folder will be deleted after test
DEBUG_MODE = False

class TestNCSDataCollection(unittest.TestCase):


  def setUp(self):
    """ Create temporary folder for temp files produced during test """
    self.current_dir = os.getcwd()
    now = datetime.now().strftime("%I%p_%m_%d_%Y")
    test_name = self.__class__.__name__
    self.tempdir = '{}_{}'.format(test_name,now)
    if not os.path.isdir(self.tempdir):
      os.mkdir(self.tempdir)
    os.chdir(self.tempdir)

  def test_look_at_file_summary(self):
    """     look at __repr__ results     """
    obj = collect_ncs_files.File_records()
    obj.data_completeness = 10
    obj.r_free_header = 0.4
    r = collect_ncs_files.Refinement_results()
    r.r_free_final  = 0.1
    r.clashscore = 10
    obj.refinement_records['init'] = r
    # print out_str

  def tearDown(self):
    """ remove temp files and folder
        When DEBUG_MODE = True, the temporary folder will not be deleted """
    if not DEBUG_MODE:
      os.chdir(self.current_dir)
      shutil.rmtree(self.tempdir)

def run_selected_tests():
  """  Run selected tests

  1) List in "tests" the names of the particular test you want to run
  2) Comment out unittest.main()
  3) Un-comment unittest.TextTestRunner().run(run_selected_tests())
  """
  tests = ['some_test_name']
  suite = unittest.TestSuite(map(MyTestCase, tests))
  return suite


if __name__ == '__main__':
  # use for individual tests
  DEBUG_MODE = True
  # unittest.TextTestRunner().run(run_selected_tests())

  # Use to run all tests
  unittest.main(verbosity=0)
