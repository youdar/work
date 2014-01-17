from __future__ import division
from iotbx import reflection_file_reader
from libtbx import easy_run
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
    self.tempdir = tempfile.mkdtemp('tempdir')
    os.chdir(self.tempdir)
    self.asu0_filename = 'asu0.pdb'
    self.ncs1_filename = 'ncs1.pdb'
    # Create and write a file ncs0.pdb
    prefix = 'ncs0'
    pdb_inp = pdb.input(source_info=None, lines=ncs0_pdb)
    pdb_inp.write_pdb_file(file_name = '{}.pdb'.format(prefix))
    xrs = pdb_inp.xray_structure_simple()
    f_obs = abs(xrs.structure_factors(d_min = 2).f_calc())
    mtz_dataset = f_obs.as_mtz_dataset(column_root_label="F-obs")
    mtz_dataset.add_miller_array(
      miller_array=f_obs.generate_r_free_flags(),
      column_root_label="R-free-flags")
    mtz_object = mtz_dataset.mtz_object()
    mtz_object.write(file_name = '{}_map.mtz'.format(prefix))
    # Process the mtz_object
    miller_arrays = mtz_object.as_miller_arrays()
    #
    #cmd = " ".join([
      #"phenix.refine",
      #"{0}.pdb {0}_map.mtz".format(prefix),
      #"strategy=none",
      #"main.number_of_macro_cycles=0",
      #"output.prefix={}".format(prefix),
      #"--overwrite",
      #"--quiet"])
    #easy_run.call(cmd)
    #miller_arrays = reflection_file_reader.any_reflection_file(
      #file_name = "{}_map.mtz".format(prefix)).as_miller_arrays()
    labels = []
    for ma in miller_arrays:
      l = ",".join(ma.info().labels)
      if(ma.is_complex_array()): labels.append(l)
    for i, l in enumerate(labels):
      cmd = " ".join([
        "phenix.map_box", pdbf, "%s_map.mtz"%prefix,
        "label='%s'"%l,
        ">zlog_%s_%s"%(str(i),prefix)])
      easy_run.call(cmd)

    print os.getcwd()


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
    ##self.assertRaises(ValueError,archive_dir.archive_dir_files,self.tempdir)
    #pass


  def tearDown(self):
    '''remove temp files and folder'''
    os.chdir(self.currnet_dir)
    shutil.rmtree(self.tempdir)


# Raw data used to build test cases
ncs0_pdb="""\n
CRYST1    9.114   10.768   11.649  90.00  90.00  90.00 P 1
SCALE1      0.109721  0.000000  0.000000        0.00000
SCALE2      0.000000  0.092868  0.000000        0.00000
SCALE3      0.000000  0.000000  0.085844        0.00000
ATOM      1  N   THR 1   1       9.483   9.256  10.995  1.00 26.11           N
ATOM      2  CA  THR 1   1       9.489   7.820  10.603  1.00 27.16           C
ATOM      3  C   THR 1   1       9.725   6.911  11.796  1.00 20.29           C
ATOM      4  O   THR 1   1      10.586   7.189  12.629  1.00 35.00           O
ATOM      5  CB  THR 1   1      10.632   7.506   9.642  1.00 34.84           C
ATOM      6  OG1 THR 1   1      10.831   8.604   8.746  1.00 67.35           O
ATOM      7  CG2 THR 1   1      10.308   6.254   8.859  1.00 43.17           C
TER
"""
asu0_pdb="""\n
CRYST1   22.845   22.121   22.495  90.00  90.00  90.00 P 1
SCALE1      0.043773  0.000000  0.000000        0.00000
SCALE2      0.000000  0.045206  0.000000        0.00000
SCALE3      0.000000  0.000000  0.044454        0.00000
ATOM      1  N   THR 1   1       9.483   9.256  10.995  1.00 26.11           N
ATOM      2  CA  THR 1   1       9.489   7.820  10.603  1.00 27.16           C
ATOM      3  C   THR 1   1       9.725   6.911  11.796  1.00 20.29           C
ATOM      4  O   THR 1   1      10.586   7.189  12.629  1.00 35.00           O
ATOM      5  CB  THR 1   1      10.632   7.506   9.642  1.00 34.84           C
ATOM      6  OG1 THR 1   1      10.831   8.604   8.746  1.00 67.35           O
ATOM      7  CG2 THR 1   1      10.308   6.254   8.859  1.00 43.17           C
TER
ATOM      1  N   THRaa   1       5.155   8.144  14.268  1.00 26.11           N
ATOM      2  CA  THRaa   1       5.854   7.737  13.019  1.00 27.16           C
ATOM      3  C   THRaa   1       7.251   7.213  13.300  1.00 20.29           C
ATOM      4  O   THRaa   1       7.985   7.794  14.097  1.00 35.00           O
ATOM      5  CB  THRaa   1       6.064   8.922  12.081  1.00 34.84           C
ATOM      6  OG1 THRaa   1       4.935   9.799  12.144  1.00 67.35           O
ATOM      7  CG2 THRaa   1       6.253   8.424  10.667  1.00 43.17           C
TER
ATOM      1  N   THRab   1       5.628   2.908  16.011  1.00 26.11           N
ATOM      2  CA  THRab   1       5.510   3.767  14.801  1.00 27.16           C
ATOM      3  C   THRab   1       6.705   4.691  14.647  1.00 20.29           C
ATOM      4  O   THRab   1       7.159   5.288  15.621  1.00 35.00           O
ATOM      5  CB  THRab   1       4.305   4.699  14.891  1.00 34.84           C
ATOM      6  OG1 THRab   1       3.216   4.028  15.533  1.00 67.35           O
ATOM      7  CG2 THRab   1       3.896   5.135  13.502  1.00 43.17           C
TER
"""


if __name__ == "__main__":
  unittest.main(verbosity=2)  # provides a command-line interface to the test script
  #unittest.main()