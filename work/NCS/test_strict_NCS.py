from __future__ import division
from iotbx.pdb.multimer_reconstruction import multimer
from iotbx import reflection_file_reader
from mmtbx.utils import flex
from libtbx import easy_run
from iotbx import pdb
from mmtbx import f_model
import unittest
import cProfile
import shutil
import tempfile
import sys,os
import time

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
    The files information was produced using
    c:\Phenix\Dev\Work\work\NCS\create_test_file.pyncs0_pdb, which contains
    a single NCS, CRYST1 and MTRIX records.
    produce the strings ncs0_pdb, asu0_pdb and ncs1_pdb for the tests below

    1) Create ncs0.pdb, with complete Asymmetric Unit (ASU)
    2) Create ncs1.pdb, pertubed version of ncs0_pdb, with a single NCS and MTRIX info
    3) Create asu0.pdb, A complete ASU produced form ncs0
    '''
    self.currnet_dir = os.getcwd()
    self.tempdir = tempfile.mkdtemp('tempdir')
    os.chdir(self.tempdir)
    self.asu0_filename = 'asu0.pdb'
    self.ncs0_filename = 'ncs0.pdb'
    self.ncs1_filename = 'ncs1.pdb'
    # Create and write a file ncs0.pdb with complete Asymmetric Unit (ASU)
    f = open(self.ncs0_filename,'w')
    f.write(ncs0_pdb)
    f.close()
    # Create and write ncs1.pdb
    f = open(self.ncs1_filename,'w')
    f.write(ncs1_pdb)
    f.close()
    # Create and write asu0.pdb
    f = open(self.asu0_filename,'w')
    f.write(asu0_pdb)
    f.close()

    # Produce experimental data for asu0
    self.pdb_inp_asu0 = pdb.input(file_name=self.asu0_filename)
    self.xrs_asu0 = self.pdb_inp_asu0.xray_structure_simple()
    self.f_obs_asu0 = abs(self.xrs_asu0.structure_factors(d_min = 2).f_calc())
    #self.f_obs_asu0 = self.xrs_asu0.structure_factors(d_min = 2).f_calc()




    #mtz_dataset = f_obs.as_mtz_dataset(column_root_label="F-obs")
    #mtz_dataset.add_miller_array(
      #miller_array=f_obs.generate_r_free_flags(),
      #column_root_label="R-free-flags")
    #mtz_object = mtz_dataset.mtz_object()
    #mtz_object.write(file_name = '{}_map.mtz'.format(prefix))
    ## Process the mtz_object
    #miller_arrays = mtz_object.as_miller_arrays()
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
    #labels = []
    #for ma in miller_arrays:
      #l = ",".join(ma.info().labels)
      #if(ma.is_complex_array()): labels.append(l)
    #for i, l in enumerate(labels):
      #cmd = " ".join([
        #"phenix.map_box", pdbf, "%s_map.mtz"%prefix,
        #"label='%s'"%l,
        #">zlog_%s_%s"%(str(i),prefix)])
      #easy_run.call(cmd)

    print os.getcwd()



  def test_process_integrity(self):
    ''' Test that when comparing asu0 to itself we get R-work is zero'''
    # Reconstruct ncs0 and retrive f_obs
    m = multimer(pdb_input_file_name=self.ncs0_filename,
                 reconstruction_type='cau',
                 error_handle=True,
                 eps=0.01)
    m.write(pdb_output_file_name='asu0_calc.pdb')
    pdb_inp = pdb.input(file_name='asu0_calc.pdb')
    xrs_calc = pdb_inp.xray_structure_simple()
    f_obs_calc = abs(xrs_calc.structure_factors(d_min = 2).f_calc())

    f1 = self.f_obs_asu0.data()
    f2 = f_obs_calc.data()

    # this way only works if f1 and f2 are of the same length
    scale = flex.sum( f1 * f2 )/ flex.sum(f2*f2)
    r_factor = flex.sum( flex.abs( f1 - scale*f2 ) ) / flex.sum( flex.abs(f1+f2) ) * 2
    self.assertEqual(r_factor,0, msg='Problem with test data, f_obs from ASU do not match those from from the same ASU as constructed by NCS')


  def test_pertubed_ncs(self):
    '''Test that the pertubed NCS (ncs1.pd) is different than the original one (ncs0_pdb)
    by checking that R-work is not zero
    Compare f_obs from asu0.pdb to f_calc from asu1.pdb'''
    # Reconstruct ncs1 and retrive f_obs
    m = multimer(pdb_input_file_name=self.ncs1_filename,
                 reconstruction_type='cau',
                 error_handle=True,
                 eps=0.01)
    m.write(pdb_output_file_name='asu1_calc.pdb')
    pdb_inp = pdb.input(file_name='asu1_calc.pdb')
    xrs_calc = pdb_inp.xray_structure_simple()
    f_obs_calc = abs(xrs_calc.structure_factors(d_min = 2).f_calc())

    f_obs_calc = f_obs_calc.common_set(self.f_obs_asu0)

    f1 = self.f_obs_asu0.data()
    f2 = f_obs_calc.data()


    scale = flex.sum( f1 * f2 )/ flex.sum(f2*f2)
    print scale

    r_factor = flex.sum( flex.abs( f1 - scale*f2 ) ) / flex.sum( flex.abs(f1+f2) ) * 2
    print r_factor
    msg='''\
    Problem with test data, f_obs from ASU do not match those from
    the same ASU as constructed by NCS'''
    self.assertTrue(r_factor>0, msg)

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

  def get_r_factor(self,f_obs,f_calc):
    # Note that f_model requires that
    # f_obs is real and that f_calc is complex
    fmodel = f_model.manager(
      f_obs = f_obs,
      f_calc = f_calc,
      f_mask = f_calc.customized_copy(data = f_calc.data()*0))
    r_factor = fmodel.r_work()
    return r_factor


  def tic(self):
    #Homemade version of matlab tic and toc functions
    global startTime_for_tictoc
    startTime_for_tictoc = time.time()

  def toc(self,msg='',print_time=True):
    if 'startTime_for_tictoc' in globals():
      if print_time:
        outstr = '{0}: Elapsed time is: {1:.4f} seconds\n'.format(msg,time.time() - startTime_for_tictoc)
        print outstr
      else:
        outstr = '{0:.4f}'.format(time.time() - startTime_for_tictoc)
        return outstr
    else:
      print "Toc: start time not set"


# Raw data used to build test cases
ncs0_pdb="""\
CRYST1   23.826   24.019   24.057  90.00  90.00  90.00 P 1
SCALE1      0.041971  0.000000  0.000000        0.00000
SCALE2      0.000000  0.041634  0.000000        0.00000
SCALE3      0.000000  0.000000  0.041568        0.00000
MTRIX1   1  1.000000  0.000000  0.000000        0.00000    1
MTRIX2   1  0.000000  1.000000  0.000000        0.00000    1
MTRIX3   1  0.000000  0.000000  1.000000        0.00000    1
MTRIX1   2  0.496590 -0.643597  0.582393        0.00000
MTRIX2   2  0.867925  0.376088 -0.324443        0.00000
MTRIX3   2 -0.010221  0.666588  0.745356        0.00000
MTRIX1   3 -0.317946 -0.173437  0.932111        0.00000
MTRIX2   3  0.760735 -0.633422  0.141629        0.00000
MTRIX3   3  0.565855  0.754120  0.333333        0.00000
ATOM      1  N   THR 1   1       9.483  10.256  10.995  1.00 26.11           N
ATOM      2  CA  THR 1   1       9.489   8.820  10.603  1.00 27.16           C
ATOM      3  C   THR 1   1       9.725   7.911  11.796  1.00 20.29           C
ATOM      4  O   THR 1   1      10.586   8.189  12.629  1.00 35.00           O
ATOM      5  CB  THR 1   1      10.632   8.506   9.642  1.00 34.84           C
ATOM      6  OG1 THR 1   1      10.831   9.604   8.746  1.00 67.35           O
ATOM      7  CG2 THR 1   1      10.308   7.254   8.859  1.00 43.17           C
"""
asu0_pdb="""\
CRYST1   23.826   24.019   24.057  90.00  90.00  90.00 P 1
SCALE1      0.041971  0.000000  0.000000        0.00000
SCALE2      0.000000  0.041634  0.000000        0.00000
SCALE3      0.000000  0.000000  0.041568        0.00000
ATOM      1  N   THR 1   1       9.483  10.256  10.995  1.00 26.11           N
ATOM      2  CA  THR 1   1       9.489   8.820  10.603  1.00 27.16           C
ATOM      3  C   THR 1   1       9.725   7.911  11.796  1.00 20.29           C
ATOM      4  O   THR 1   1      10.586   8.189  12.629  1.00 35.00           O
ATOM      5  CB  THR 1   1      10.632   8.506   9.642  1.00 34.84           C
ATOM      6  OG1 THR 1   1      10.831   9.604   8.746  1.00 67.35           O
ATOM      7  CG2 THR 1   1      10.308   7.254   8.859  1.00 43.17           C
TER
ATOM      1  N   THRaa   1       4.512   8.520  14.935  1.00 26.11           N
ATOM      2  CA  THRaa   1       5.211   8.113  13.685  1.00 27.16           C
ATOM      3  C   THRaa   1       6.608   7.589  13.966  1.00 20.29           C
ATOM      4  O   THRaa   1       7.342   8.170  14.764  1.00 35.00           O
ATOM      5  CB  THRaa   1       5.421   9.299  12.748  1.00 34.84           C
ATOM      6  OG1 THRaa   1       4.291  10.175  12.810  1.00 67.35           O
ATOM      7  CG2 THRaa   1       5.610   8.800  11.333  1.00 43.17           C
TER
ATOM      1  N   THRab   1       5.455   2.275  16.765  1.00 26.11           N
ATOM      2  CA  THRab   1       5.336   3.134  15.555  1.00 27.16           C
ATOM      3  C   THRab   1       6.531   4.058  15.401  1.00 20.29           C
ATOM      4  O   THRab   1       6.986   4.655  16.375  1.00 35.00           O
ATOM      5  CB  THRab   1       4.132   4.066  15.645  1.00 34.84           C
ATOM      6  OG1 THRab   1       3.043   3.395  16.287  1.00 67.35           O
ATOM      7  CG2 THRab   1       3.722   4.502  14.256  1.00 43.17           C
TER
"""

ncs1_pdb = """\
CRYST1   23.826   24.019   24.057  90.00  90.00  90.00 P 1
SCALE1      0.041971  0.000000  0.000000        0.00000
SCALE2      0.000000  0.041634  0.000000        0.00000
SCALE3      0.000000  0.000000  0.041568        0.00000
MTRIX1   1  1.000000  0.000000  0.000000        0.00000    1
MTRIX2   1  0.000000  1.000000  0.000000        0.00000    1
MTRIX3   1  0.000000  0.000000  1.000000        0.00000    1
MTRIX1   2  0.496590 -0.643597  0.582393        0.00000
MTRIX2   2  0.867925  0.376088 -0.324443        0.00000
MTRIX3   2 -0.010221  0.666588  0.745356        0.00000
MTRIX1   3 -0.317946 -0.173437  0.932111        0.00000
MTRIX2   3  0.760735 -0.633422  0.141629        0.00000
MTRIX3   3  0.565855  0.754120  0.333333        0.00000
ATOM      1  N   THR 1   1       9.514  10.278  10.935  1.00 26.11           N
ATOM      2  CA  THR 1   1       9.519   8.854  10.711  1.00 27.16           C
ATOM      3  C   THR 1   1       9.699   7.867  11.829  1.00 20.29           C
ATOM      4  O   THR 1   1      10.692   8.121  12.544  1.00 35.00           O
ATOM      5  CB  THR 1   1      10.676   8.475   9.646  1.00 34.84           C
ATOM      6  OG1 THR 1   1      10.859   9.626   8.815  1.00 67.35           O
ATOM      7  CG2 THR 1   1      10.262   7.276   8.864  1.00 43.17           C
TER
"""

if __name__ == "__main__":
  unittest.main(verbosity=2)  # provides a command-line interface to the test script
  #unittest.main()