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
    create a temporary folder with files for testing, based
    on ncs_0_copy (see below)

    1) Create ncs0.pdb, with complete Asymmetric Unit (ASU)
    2) Create ncs1.pdb, pertubed (shaken) version of ncs0_pdb,
       with a single NCS and MTRIX info
    3) Create asu0.pdb, A complete ASU produced form ncs0
    4) Create asu1.pdb, A complete ASU produced form ncs1
    '''
    self.currnet_dir = os.getcwd()
    self.tempdir = tempfile.mkdtemp('tempdir')
    # for testing use junk folder insted of real temp directory

    # remember to change back the clean up when going back to real temp dir
    osType = sys.platform
    if osType.startswith('win'):
      self.tempdir = (r'C:\Phenix\Dev\Work\work\NCS\junk')
    else:
      self.tempdir = ('/net/cci/youval/Work/work/NCS/junk')

    os.chdir(self.tempdir)
    ncs0_filename = 'ncs0.pdb'
    ncs1_filename = 'ncs1_shaken.pdb'
    asu0_filename = 'asu0.pdb'
    asu1_filename = 'asu1_shaken.pdb'

    # Write single copy NCS file with MTRIX record
    open(ncs0_filename,'w').write(ncs_0_copy)

    # Create the ASU coordinates using MTRIX records
    crystal_symmetry = self.create_asu(
      ncs_filename=ncs0_filename,
      asu_filename=asu0_filename)

    # Create and write a file ncs0.pdb with complete Asymmetric Unit (ASU)
    # 1 NCS copy: starting template to generate whole asu; place into P1 box
    pdb_inp_ncs = pdb.input(source_info=None, lines=ncs_0_copy)
    pdb_inp_ncs.write_pdb_file(
      file_name=ncs0_filename,
      crystal_symmetry = crystal_symmetry)

    # When using pdb_inp.write_pdb_file the MTRIX record are omitted. Add them back
    self.add_MTRIX_to_pdb(ncs0_filename, record=ncs_0_copy)

    # Shake structure - subject to refinement input
    pdb_inp_ncs = pdb.input(file_name = ncs0_filename)
    # Create xrs from the single ncs
    xrs = pdb_inp_ncs.xray_structure_simple()
    xrs_one_ncs = xrs.orthorhombic_unit_cell_around_centered_scatterers(
      buffer_size=8)
    ph = pdb_inp_ncs.construct_hierarchy()
    ph.adopt_xray_structure(xrs)
    xrs_shaken = xrs_one_ncs.deep_copy_scatterers()
    xrs_shaken.shake_sites_in_place(mean_distance=0.3)
    ph.adopt_xray_structure(xrs_shaken)
    ph.write_pdb_file(file_name=ncs1_filename, crystal_symmetry = crystal_symmetry)
    self.add_MTRIX_to_pdb(ncs0_filename, record=ncs_0_copy)


    m = multimer(ncs_filename,'cau',error_handle=True,eps=1e-2)
    m.write(self.asu0_filename)
    pdb_inp_asu = pdb.input(file_name = asu_filename)
    # Create xrs from the complete asu
    xrs = pdb_inp_asu.xray_structure_simple()
    crystal_symmetry = xrs.crystal_symmetry()
    # write the pdb files again with MTRIX and CRYST1 records
    mtrix_object = pdb_inp_ncs.process_mtrix_records()

    pdb_inp.write_pdb_file(file_name = asu_filename, crystal_symmetry = crystal_symmetry)

    mtrix_object = pdb_inp_ncs.process_mtrix_records()

    # 1 NCS copy -> full asu (expand NCS). This is the answer-strucure
    #m = multimer("one_ncs_in_asu.pdb",'cau',error_handle=True,eps=1e-2)
    #assert m.number_of_transforms == 3

    ph = pdb_inp.construct_hierarchy()
    xrs = pdb_inp.xray_structure_simple()
    xrs_one_ncs = xrs.orthorhombic_unit_cell_around_centered_scatterers(
      buffer_size=8)
    ph.adopt_xray_structure(xrs)
    of = open("one_ncs_in_asu.pdb", "w")
    print >> of, mtrix_object.format_MTRIX_pdb_string()
    print >> of, ph.as_pdb_string(crystal_symmetry=xrs.crystal_symmetry())
    of.close()


    #f = open(self.ncs0_filename,'w').write(ncs0_pdb)
    ## Create and write ncs1.pdb
    #f = open(self.ncs1_filename,'w').write(ncs1_pdb)
    ## Create and write asu0.pdb
    #f = open(self.asu0_filename,'w').write(asu0_pdb)

    ## Produce experimental data for asu0
    #self.pdb_inp_asu0 = pdb.input(file_name=self.asu0_filename)
    #self.xrs_asu0 = self.pdb_inp_asu0.xray_structure_simple()
    #self.f_obs_asu0 = abs(self.xrs_asu0.structure_factors(d_min = 2).f_calc())
    ##

    pdb_inp = pdb.input(source_info=None, lines=ncs_0_copy)
    mtrix_object = pdb_inp.process_mtrix_records()
    ph = pdb_inp.construct_hierarchy()
    xrs = pdb_inp.xray_structure_simple()
    xrs_one_ncs = xrs.orthorhombic_unit_cell_around_centered_scatterers(
      buffer_size=8)
    ph.adopt_xray_structure(xrs)
    of = open("one_ncs_in_asu.pdb", "w")
    print >> of, mtrix_object.format_MTRIX_pdb_string()
    print >> of, ph.as_pdb_string(crystal_symmetry=xrs.crystal_symmetry())
    of.close()


    # 1 NCS copy -> full asu (expand NCS). This is the answer-strucure
    m = multimer("one_ncs_in_asu.pdb",'cau',error_handle=True,eps=1e-2)
    assert m.number_of_transforms == 3
    xrs_asu = m.assembled_multimer.extract_xray_structure(
      crystal_symmetry = xrs.crystal_symmetry())
    m.write("full_asu.pdb")
    # Shake structure - subject to refinement input
    xrs_shaken = xrs_one_ncs.deep_copy_scatterers()
    xrs_shaken.shake_sites_in_place(mean_distance=0.3)
    ph.adopt_xray_structure(xrs_shaken)
    of = open("one_ncs_in_asu_shaken.pdb", "w")
    print >> of, mtrix_object.format_MTRIX_pdb_string()
    print >> of, ph.as_pdb_string(crystal_symmetry=xrs.crystal_symmetry())



    # Generate Fobs from answer structure
    f_obs = abs(xrs_asu.structure_factors(d_min=2, algorithm="direct").f_calc())
    r_free_flags = f_obs.generate_r_free_flags()
    mtz_dataset = f_obs.as_mtz_dataset(column_root_label="F-obs")
    mtz_dataset.add_miller_array(
      miller_array=r_free_flags,
      column_root_label="R-free-flags")
    mtz_object = mtz_dataset.mtz_object()
    mtz_object.write(file_name = "data.mtz")

    # Refinement


  def test_process_integrity(self):
    ''' Test that when comparing asu0 to itself we get R-work is zero'''
    # Reconstruct ncs0 and retrive f_obs
    m = self.build_asu(file_name_ncs=self.ncs0_filename, file_name_asu='asu0.pdb')
    f_calc = self.get_f_calc(file_name='asu0.pdb')
    r_factor = self.calc_r_factor(self.f_obs_asu0,f_calc)
    msg='Problem with test data, f_obs from ASU do not match those from from the same ASU as constructed by NCS'
    self.assertEqual(r_factor,0, msg)

  def test_pertubed_ncs(self):
    '''
    Test that the perturbed  NCS (ncs1.pdb) is different than the original one (ncs0_pdb)
    by checking that R-work > 0.3
    Compare f_obs from asu0.pdb to f_calc from asu1.pdb'''
    # Reconstruct ncs1 and retrive f_obs
    m = self.build_asu(file_name_ncs=self.ncs1_filename, file_name_asu='asu1.pdb')
    f_calc = self.get_f_calc(file_name='asu1.pdb')
    r_factor = self.calc_r_factor(self.f_obs_asu0,f_calc)
    msg='''\
    Problem with test data, f_obs from ASU do not match those from
    the same ASU as constructed by NCS'''
    self.assertTrue(r_factor > 0.3, msg)

  def test_refinement(self):
    '''
    Test that refining asu1.pdb, build from the perturbed NCS, converge to asu0.pdb
    Use asu0.pdb to create x-ray structure (xrs)
    and from it produce F_obs (observed structure factors)'''
    # produce the ASUs from the NCSs, in the current directory
    m_asu0 = self.build_asu(file_name_ncs=self.ncs0_filename, file_name_asu=self.asu0_filename)
    m_asu1 = self.build_asu(file_name_ncs=self.ncs1_filename, file_name_asu=self.asu1_filename)
    #
    f_obs = abs(self.get_f_calc(file_name=self.asu0_filename))
    f_calc = self.get_f_calc(file_name=self.asu1_filename)
    # r_factor at start
    r_factor = self.calc_r_factor(f_obs,f_calc)
    # Refine
    # 1) create mtz file
    file_name_mtz = self.asu0_filename.split('.')[0] + '_map.mtz'
    self.make_mtz_file(f_obs=f_obs, file_name_mtz=file_name_mtz)
    # 2) Refine
    self.call_refine(
      pdb_file=self.asu1_filename,
      mtz_file=file_name_mtz,
      output_file_name='refine_output')
      #pdb_file_symmetry_target=self.asu0_filename)
    # Process refinement resaults
    file_name_refined = 'refine_output_001.pdb'
    self.assertTrue(os.path.isfile(file_name_refined),'No {} refined pdb file'.format(file_name_refined))
    f_calc = self.get_f_calc(file_name=file_name_refined)
    r_factor_refined = self.get_r_factor(f_obs=f_obs,f_calc=f_calc)
    msg1='r_factor is {0:.5f}. very small before refinement.'.format(r_factor)
    self.assertTrue(r_factor>0.3,msg1)
    msg2='Refinement did not work well, r_factor before {0:.5f}  after {1:.5f}'.format(r_factor,r_factor_refined)
    self.assertTrue(r_factor_refined<0.05,msg2)

  def test_ncs_refinement(self):
    '''Test refinement using strict NCS.
    Refinement using the gradient of only one NCS copy'''
    pass

  #def test_2(self):
    #'''Test that Error is raised when trying to ...'''
    ##self.assertRaises(ValueError,..)
    #pass

  def tearDown(self):
    '''remove temp files and folder'''
    #os.chdir(self.currnet_dir)
    #shutil.rmtree(self.tempdir)

  def create_asu(self,ncs_filename,asu_filename,crystal_symmetry=None):
    ''' (str,str) -> crystal_symmetry object
    Create ASU from NCS and save the new pdb file
    with CRYST1 and SCALE records in local directory

    Argument:
    ---------
    ncs_filename : NCS to be used to create the ASU
    asu_filename : ASU file name
    crystal_symmetry : Allow forcing other crystal_symmetry

    Returns:
    --------
    crystal_symmetry
    '''
    m = multimer(ncs_filename,'cau',error_handle=True,eps=1e-2)
    if m.number_of_transforms == 0:
      print 'Number of transforms is zero'
    m.write(asu_filename)
    pdb_inp = pdb.input(file_name = asu_filename)
    xrs = pdb_inp.xray_structure_simple()
    if not crystal_symmetry:
      crystal_symmetry = xrs.crystal_symmetry()
    pdb_inp.write_pdb_file(file_name = asu_filename, crystal_symmetry = crystal_symmetry)
    return crystal_symmetry

  def add_MTRIX_to_pdb(self,pdb_fn, record):
    '''(str,ste) -> None
    Add MTRIX records from the string record
    to the file pdb_fn and write the modified pdb file, with
    the MTRIX records, in the current directory
    '''
    ncs_pdb = open(pdb_fn,'r').read().splitlines()
    record = record.splitlines()
    i = 0
    while record:
      x = record.pop(0)
      if x.startswith('MTRIX'):
        ncs_pdb.insert(i+4, x)
        i += 1
    ncs_pdb = '\n'.join(ncs_pdb)
    open(pdb_fn,'w').write(ncs_pdb)

  def call_refine(self,pdb_file,mtz_file,output_file_name,pdb_file_symmetry_target=None):
    '''
    Run refinement and produce refined pdb file in current directory

    To set the number of refinement cycles change main.number_of_macro_cycles

    Argument:
    ---------
    pdb_file : pdb file to be refined
    mtz_file : mtz file from the target experiment f_obs
    output_file_name : the output file prefix
    pdb_file_symmetry_target : use if you want to force symmetry of a pdb file

    Output:
    -------
    writing out refined model, complete refinement statistics and
    electron density maps in various formats.
    '''
    cmd = " ".join([
      "phenix.refine",
      "{0} {1}".format(pdb_file,mtz_file),
      "strategy=individual_sites",
      "main.number_of_macro_cycles=3",
      "output.prefix={}".format(output_file_name),
      "--overwrite",
      "--quiet"])
    if pdb_file_symmetry_target:
      cmd = ' '.join([cmd,"--symmetry={}".format(pdb_file_symmetry_target)])
    easy_run.call(cmd)

  def make_mtz_file(self,f_obs,file_name_mtz):
    '''
    Create mtz file, from f_obs or f_calc, in the current directory
    file_name_mtz need to have the format name.mtz
    '''
    mtz_dataset = abs(f_obs).as_mtz_dataset(column_root_label="F-obs")
    mtz_dataset.add_miller_array(
      miller_array=f_obs.generate_r_free_flags(),
      column_root_label="R-free-flags")
    mtz_object = mtz_dataset.mtz_object()
    mtz_object.write(file_name=file_name_mtz)

  def get_f_calc(self,file_name):
    '''Calculate f_calc (diffraction image frequencies) from a pdb file'''
    pdb_inp = pdb.input(file_name=file_name)
    xrs = pdb_inp.xray_structure_simple()
    f_calc = xrs.structure_factors(d_min = 2).f_calc()
    return f_calc

  def build_asu(self,file_name_ncs,file_name_asu):
    '''Build ASU from NCS and save the new pdb file in local directory'''
    m = multimer(
      pdb_input_file_name=file_name_ncs,
      reconstruction_type='cau',
      error_handle=True,
      eps=0.01)
    m.write(pdb_output_file_name=file_name_asu)
    return m

  def calc_r_factor(self,f_obs,f_calc):
    ''' Both f_obs and f_calc need to be real'''
    f1 = abs(f_obs).data()
    f2 = abs(f_calc).data()
    scale = flex.sum( f1 * f2 )/ flex.sum(f2*f2)
    r_factor = flex.sum( flex.abs( f1 - scale*f2 ) ) / flex.sum( flex.abs(f1+f2) ) * 2
    return r_factor

  def get_r_factor(self,f_obs,f_calc):
    '''
    When using f_model to get r_factor,
    Note that f_model requires that
    f_obs is real and that f_calc is complex
    '''
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
ncs_0_copy="""\
MTRIX1   1  1.000000  0.000000  0.000000        0.00000    1
MTRIX2   1  0.000000  1.000000  0.000000        0.00000    1
MTRIX3   1  0.000000  0.000000  1.000000        0.00000    1
MTRIX1   2  0.496590 -0.643597  0.582393        0.00000    1
MTRIX2   2  0.867925  0.376088 -0.324443        0.00000    1
MTRIX3   2 -0.010221  0.666588  0.745356        0.00000    1
MTRIX1   3 -0.317946 -0.173437  0.932111        0.00000
MTRIX2   3  0.760735 -0.633422  0.141629        0.00000
MTRIX3   3  0.565855  0.754120  0.333333        0.00000
ATOM      1  N   THR A   1       9.670  10.289  11.135  1.00 20.00           N
ATOM      2  CA  THR A   1       9.559   8.931  10.615  1.00 20.00           C
ATOM      3  C   THR A   1       9.634   7.903  11.739  1.00 20.00           C
ATOM      4  O   THR A   1      10.449   8.027  12.653  1.00 20.00           O
ATOM      5  CB  THR A   1      10.660   8.630   9.582  1.00 20.00           C
ATOM      6  OG1 THR A   1      10.560   9.552   8.490  1.00 20.00           O
ATOM      7  CG2 THR A   1      10.523   7.209   9.055  1.00 20.00           C
TER
"""

if __name__ == "__main__":
  unittest.main(verbosity=2)  # provides a command-line interface to the test script
  #unittest.main()


'''
    ## Process the mtz_object
    #miller_arrays = mtz_object.as_miller_arrays()
    #
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
'''
