from __future__ import division
from datetime import datetime
import collect_ncs_files
import cPickle as pickle
import unittest
import get_mtz
import shutil
import sys
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
    print obj

  def test_getting_file_info(self):
    """ making sure MTRIX records are handled correctly   """
    pdbs = [test_pdb_1,test_pdb_2,test_pdb_3,test_pdb_4]
    fns = ['test_pdb_1','test_pdb_2','test_pdb_3','test_pdb_4']
    n_ncs = [2,None,3,None]
    for fn,s,n in zip(fns,pdbs,n_ncs):
      open(fn+'.pdb','w').write(s)
      collect = collect_ncs_files.ncs_paper_data_collection()
      a = collect.get_pdb_file_info(fn)
      if a:
        n_ncs_copies = a.n_ncs_copies
      else:
        n_ncs_copies = None
      self.assertEqual(n_ncs_copies,n)

    # clean test pdb files from pdb folder
    files_to_delete = [
      'test_pdb_1.pdb','test_pdb_2.pdb','test_pdb_3.pdb','test_pdb_4.pdb']
    for fn in files_to_delete:
      full_path = os.path.join(collect.asu_dir,fn)
      if os.path.isfile(full_path):
        os.remove(full_path)

  def test_reading_processing_pdb_file(self):
    """ test processing a pdb file """
    collect = collect_ncs_files.ncs_paper_data_collection()
    a = collect.get_pdb_file_info('1vcr')
    self.assertEqual(a.n_ncs_copies,5)

  def test_command_line(self):
    """     Test the command line call    """
    osType = sys.platform
    if osType.startswith('win'):
      pass
    else:
      sources = os.environ["workfolder"]
      com_path = sources + '/NCS/ncs_paper/get_file_ncs_info.py'
      cmd = 'python {} 1vcr > log_test'.format(com_path)
      os.system(cmd)

  def test_reading_results(self):
    """ make sure we can read the results """
    sources = os.environ["workfolder"]
    path = sources + '/NCS/ncs_paper/ncs_queue_results'
    fn = os.path.join(path,'log_1vcr')
    if os.path.isfile(fn):
      pdb_info = pickle.load(open(fn,'r'))
      print pdb_info
    else:
      print 'Could not find log_vcr'

  def test_records_update(self):
    """ test that we reformat the file records object
    attribute names correctly """
    class Old_File_records(object):

      def __init__(self):

        self.pdb_id = ''
        self.n_ncs_copies = None
        self.year = None
        self.resolution = None
        self.data_completeness = None
        self.solvent = None
        self.only_master_in_pdb = None
        self.ncs_reported_in_pdb = None
        # data_to_param_ratio: f_obs.size / atom_number
        self.data_to_param_ration = None
        # refinement_records contains Refinement_results objects
        # for example refinement_records['using cartesian NCS']
        self.refinement_records = {}

    f = Old_File_records()
    f.solvent = 7
    f.n_ncs_copies = 10
    f.refinement_records['xray'] = [1,2,3]
    n = get_mtz.fix_pdb_info(f)
    self.assertEqual(n.solvent_fraction,7)
    self.assertEqual(n.n_ncs_copies,10)
    self.assertEqual(n.refinement_records['xray'],[1,2,3])

  def test_get_mtz(self):
    """ make sure get_mtz works """
    get_mtz.run(['1vcr'])


  def tearDown(self):
    """ remove temp files and folder
        When DEBUG_MODE = True, the temporary folder will not be deleted """
    os.chdir(self.current_dir)
    try:
      if not DEBUG_MODE:
        shutil.rmtree(self.tempdir)
    except:
      pass

test_pdb_1 = '''\
CRYST1  577.812  448.715  468.790  90.00  90.00  90.00 P 1
SCALE1      0.001731  0.000000  0.000000        0.00000
SCALE2      0.000000  0.002229  0.000000        0.00000
SCALE3      0.000000  0.000000  0.002133        0.00000
ATOM      1  CA  LYS A 151      10.766   9.333  12.905  1.00 44.22           C
ATOM      2  CA  LYS A 152      10.117   9.159  11.610  1.00 49.42           C
ATOM      3  CA  LYS A 153       9.099   8.000  11.562  1.00 46.15           C
ATOM      4  CA  LYS A 154       8.000   8.202  11.065  1.00 52.97           C
ATOM      5  CA  LYS A 155      11.146   9.065  10.474  1.00 41.68           C
ATOM      6  CA  LYS A 156      10.547   9.007   9.084  1.00 55.55           C
ATOM      7  CA  LYS A 157      11.545   9.413   8.000  1.00 72.27           C
ATOM      8  CA  LYS A 158      12.277  10.718   8.343  1.00 75.78           C
ATOM      9  CA  LYS A 159      11.349  11.791   8.809  1.00 75.88           C
TER
ATOM    222  CA  LEU X  40      94.618  -5.253  91.582  1.00 87.10           C
ATOM    223  CA  ARG X  41      62.395  51.344  80.786  1.00107.25           C
ATOM    224  CA  ARG X  42      62.395  41.344  80.786  1.00107.25           C
TER
ATOM      1  CA  THR D   1       8.111  11.080  10.645  1.00 20.00           C
ATOM      2  CA  THR D   2       8.000   9.722  10.125  1.00 20.00           C
ATOM      3  CA  THR D   3       8.075   8.694  11.249  1.00 20.00           C
ATOM      4  CA  THR D   4       8.890   8.818  12.163  1.00 20.00           C
ATOM      5  CA  THR D   5       9.101   9.421   9.092  1.00 20.00           C
ATOM      6  CA  THR D   6       9.001  10.343   8.000  1.00 20.00           C
ATOM      7  CA  THR D   7       8.964   8.000   8.565  1.00 20.00           C
TER
ATOM      1  CA  LYS B 151       6.855   8.667  15.730  1.00 44.22           C
ATOM      2  CA  LYS B 152       5.891   8.459  14.655  1.00 49.42           C
ATOM      3  CA  LYS B 153       6.103   7.155  13.858  1.00 46.15           C
ATOM      4  CA  LYS B 154       5.138   6.438  13.633  1.00 52.97           C
ATOM      5  CA  LYS B 155       5.801   9.685  13.736  1.00 41.68           C
ATOM      6  CA  LYS B 156       4.731   9.594  12.667  1.00 55.55           C
ATOM      7  CA  LYS B 157       4.334  10.965  12.119  1.00 72.27           C
ATOM      8  CA  LYS B 158       4.057  11.980  13.238  1.00 75.78           C
ATOM      9  CA  LYS B 159       3.177  11.427  14.310  1.00 75.88           C
END
'''

test_pdb_2 = '''\
CRYST1  577.812  448.715  468.790  90.00  90.00  90.00 P 1
SCALE1      0.001731  0.000000  0.000000        0.00000
SCALE2      0.000000  0.002229  0.000000        0.00000
SCALE3      0.000000  0.000000  0.002133        0.00000
ATOM      1  CA  LYS A 151      10.766   9.333  12.905  1.00 44.22           C
ATOM      2  CA  LYS A 152      10.117   9.159  11.610  1.00 49.42           C
ATOM      3  CA  LYS A 153       9.099   8.000  11.562  1.00 46.15           C
ATOM      4  CA  LYS A 154       8.000   8.202  11.065  1.00 52.97           C
ATOM      5  CA  LYS A 155      11.146   9.065  10.474  1.00 41.68           C
ATOM      6  CA  LYS A 156      10.547   9.007   9.084  1.00 55.55           C
ATOM      7  CA  LYS A 157      11.545   9.413   8.000  1.00 72.27           C
ATOM      8  CA  LYS A 158      12.277  10.718   8.343  1.00 75.78           C
ATOM      9  CA  LYS A 159      11.349  11.791   8.809  1.00 75.88           C
TER
ATOM    222  CA  LEU X  40      94.618  -5.253  91.582  1.00 87.10           C
ATOM    223  CA  ARG X  41      62.395  51.344  80.786  1.00107.25           C
ATOM    224  CA  ARG X  42      62.395  41.344  80.786  1.00107.25           C
END
'''

test_pdb_3 = """\
CRYST1  577.812  448.715  468.790  90.00  90.00  90.00 P 1
MTRIX1   1  1.000000  0.000000  0.000000        0.00000    1
MTRIX2   1  0.000000  1.000000  0.000000        0.00000    1
MTRIX3   1  0.000000  0.000000  1.000000        0.00000    1
MTRIX1   2  0.496590 -0.643597  0.582393        0.00000
MTRIX2   2  0.867925  0.376088 -0.324443        0.00000
MTRIX3   2 -0.010221  0.666588  0.745356        0.00000
MTRIX1   3 -0.317946 -0.173437  0.932111        0.00000
MTRIX2   3  0.760735 -0.633422  0.141629        0.00000
MTRIX3   3  0.565855  0.754120  0.333333        0.00000
ATOM    757  N   ASP A 247      16.068 -20.882 -28.984  1.00 35.93           N
ATOM    758  CA  ASP A 247      15.914 -22.265 -28.600  1.00 47.90           C
ATOM    759  C   ASP A 247      17.130 -23.042 -29.116  1.00 42.32           C
ATOM    760  O   ASP A 247      17.461 -22.986 -30.301  1.00 47.25           O
ATOM    761  CB  ASP A 247      14.621 -22.814 -29.198  1.00 47.22           C
ATOM    762  CG  ASP A 247      14.068 -23.974 -28.412  1.00 61.15           C
ATOM    763  OD1 ASP A 247      14.359 -24.061 -27.196  1.00 63.66           O
ATOM    764  OD2 ASP A 247      13.341 -24.798 -29.012  1.00 77.01           O
ATOM    765  N   VAL A 248      17.808 -23.746 -28.218  1.00 44.08           N
ATOM    766  CA  VAL A 248      19.008 -24.503 -28.584  1.00 46.18           C
ATOM    767  C   VAL A 248      18.668 -25.988 -28.583  1.00 53.97           C
ATOM    768  O   VAL A 248      18.049 -26.478 -27.638  1.00 51.48           O
ATOM    769  CB  VAL A 248      20.185 -24.226 -27.608  1.00 47.55           C
ATOM    770  CG1 VAL A 248      21.414 -25.015 -28.012  1.00 41.43           C
ATOM    771  CG2 VAL A 248      20.513 -22.743 -27.567  1.00 41.64           C
ATOM    772  N   VAL A 249      19.057 -26.697 -29.641  1.00 54.29           N
ATOM    773  CA  VAL A 249      18.662 -28.097 -29.810  1.00 60.17           C
ATOM    774  C   VAL A 249      19.859 -29.041 -29.982  1.00 57.98           C
ATOM    775  O   VAL A 249      20.731 -28.827 -30.828  1.00 58.31           O
ATOM    776  CB  VAL A 249      17.671 -28.280 -30.997  1.00 60.85           C
ATOM    777  CG1 VAL A 249      16.500 -27.300 -30.884  1.00 48.00           C
ATOM    778  CG2 VAL A 249      18.386 -28.110 -32.337  1.00 59.99           C
TER
ATOM    780  N   LYS D 151       4.045  -6.858 -32.823  1.00 45.22           N
ATOM    781  CA  LYS D 151       4.686  -6.715 -34.123  1.00 50.40           C
ATOM    782  C   LYS D 151       5.707  -5.554 -34.172  1.00 47.13           C
ATOM    783  O   LYS D 151       6.820  -5.764 -34.625  1.00 52.91           O
ATOM    784  CB  LYS D 151       3.657  -6.646 -35.268  1.00 40.73           C
ATOM    785  CG  LYS D 151       4.264  -6.627 -36.661  1.00 55.98           C
ATOM    786  CD  LYS D 151       3.272  -7.051 -37.745  1.00 72.14           C
ATOM    787  CE  LYS D 151       2.529  -8.338 -37.375  1.00 75.11           C
ATOM    788  NZ  LYS D 151       3.451  -9.400 -36.884  1.00 75.46           N
ATOM    789  N   ARG D 152       5.369  -4.349 -33.709  1.00 42.01           N
ATOM    790  CA  ARG D 152       6.399  -3.290 -33.702  1.00 40.51           C
ATOM    791  C   ARG D 152       6.155  -2.002 -32.909  1.00 34.21           C
ATOM    792  O   ARG D 152       5.015  -1.605 -32.636  1.00 33.77           O
END
"""

test_pdb_4 = """\
MTRIX1   1  1.000000  0.000000  0.000000        0.00000    1
MTRIX2   1  0.000000  1.000000  0.000000        0.00000    1
MTRIX3   1  0.000000  0.000000  1.000000        0.00000    1
MTRIX1   2  0.496590 -0.643597  0.582393        0.00000    1
MTRIX2   2  0.867925  0.376088 -0.324443        0.00000    1
MTRIX3   2 -0.010221  0.666588  0.745356        0.00000    1
MTRIX1   3 -0.317946 -0.173437  0.932111        0.00000    1
MTRIX2   3  0.760735 -0.633422  0.141629        0.00000    1
MTRIX3   3  0.565855  0.754120  0.333333        0.00000    1
ATOM      1  N   THR A   1       9.670  10.289  11.135  1.00 20.00           N
ATOM      2  CA  THR A   2       9.559   8.931  10.615  1.00 20.00           C
ATOM      3  C   THR A   3       9.634   7.903  11.739  1.00 20.00           C
ATOM      4  O   THR B   4      10.449   8.027  12.653  1.00 20.00           O
ATOM      5  CB  THR B   5      10.660   8.630   9.582  1.00 20.00           C
ATOM      6  OG1 THR A   6      10.560   9.552   8.490  1.00 20.00           O
ATOM      7  CG2 THR A   7      10.523   7.209   9.055  1.00 20.00           C
END
"""

def run_selected_tests():
  """  Run selected tests

  1) List in "tests" the names of the particular test you want to run
  2) Comment out unittest.main()
  3) Un-comment unittest.TextTestRunner().run(run_selected_tests())
  """
  tests = ['test_get_mtz']
  suite = unittest.TestSuite(map(TestNCSDataCollection, tests))
  return suite


if __name__ == '__main__':
  # use for individual tests
  # DEBUG_MODE = True
  unittest.TextTestRunner().run(run_selected_tests())

  # Use to run all tests
  # unittest.main(verbosity=0)
