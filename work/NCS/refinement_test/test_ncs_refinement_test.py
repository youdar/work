from __future__ import division
# from test_ncs_refinement import ncs_refine_test
from misc_scripts.strict_ncs.test_ncs_refinement import ncs_refine_test
import os
import sys


__author__ = 'Youval'

"""
This is a test program for test_ncs_refinement.py
"""

class test_ncs_refinement(object):

  def setUp(self):
    self.currnet_dir = os.getcwd()
    self.path =os.path.dirname(__file__)
    print self.path
    os.chdir(self.path)

  def test_sites(self):
    """
    test sites refinement, without geometry restraints

    Uses the files xxxx.pdb and xxxx-sf.cif, where site were shaken
    """
    print 'Running ',sys._getframe().f_code.co_name
    test_obj = ncs_refine_test(
      use_geometry_restraints=False, real_data=False,
      sf_and_grads_algorithm = 'direct',
      target_name = 'ls_wunit_k1')
    test_obj.set_working_path(self.path)
    test_obj.process_pdb_and_cif_files('xxxx')
    test_obj.refinement_loop(
    n_macro_cycle        = 100,
    r_work_target        = 0.00001,
    print_during_refinement = True,
    finite_grad_differences_test = True,
    sites                = True,
    u_iso                = False,
    transformations      = False)
    print test_obj
    print '='*173


  def test_adp(self):
    """
    test b-factor (ADP) refinement, without geometry restraints

    Uses the files yyyy.pdb and yyyy-sf.cif, where b-factors were shaken
    """
    print 'Running ',sys._getframe().f_code.co_name
    test_obj = ncs_refine_test(
      use_geometry_restraints=False, real_data=False,
      sf_and_grads_algorithm = 'direct',
      target_name = 'ls_wunit_k1')
    test_obj.set_working_path(self.path)
    test_obj.process_pdb_and_cif_files('yyyy')
    test_obj.refinement_loop(
    n_macro_cycle        = 50,
    r_work_target        = 0.00001,
    print_during_refinement = True,
    sites                = False,
    u_iso                = True,
    transformations      = False)
    print test_obj
    print '='*173

  def test_alternate(self):
    """
    test alternate sites / ADP refinement, without geometry restraints

    Uses the files zzzz.pdb and zzzz-sf.cif, where both
    site and b-factors were shaken
    """
    print 'Running ',sys._getframe().f_code.co_name
    test_obj = ncs_refine_test(
      use_geometry_restraints=False, real_data=False,
      sf_and_grads_algorithm = 'direct',
      target_name = 'ls_wunit_k1')
    test_obj.set_working_path(self.path)
    test_obj.process_pdb_and_cif_files('zzzz')
    test_obj.refinement_loop(
    n_macro_cycle        = 50,
    r_work_target        = 0.00001,
    print_during_refinement = True,
    sites                = True,
    u_iso                = True,
    transformations      = True)
    print test_obj
    print '='*173

  def test_with_geometry_restraints(self):
    """
    test sites refinement, with geometry restraints

    Uses the files xxxx.pdb and xxxx-sf.cif, where site were shaken
    """
    print 'Running ',sys._getframe().f_code.co_name
    test_obj = ncs_refine_test(
      use_geometry_restraints=True,
      real_data=False,
      sf_and_grads_algorithm = 'direct',
      target_name = 'ls_wunit_k1')
    test_obj.set_working_path(self.path)
    test_obj.process_pdb_and_cif_files('xxxx')
    test_obj.refinement_loop(
    n_macro_cycle        = 50,
    r_work_target        = 0.0001,
    print_during_refinement = True,
    sites                = True,
    u_iso                = False,
    transformations      = False)
    print test_obj
    print '='*173

  def test_transformation_refinement(self):
    """
    Test transformation refinement
    """
    print 'Running ',sys._getframe().f_code.co_name
    test_obj = ncs_refine_test(
      use_geometry_restraints=False,  # change to True
      real_data=False,
      sf_and_grads_algorithm = 'direct',
      target_name = 'ls_wunit_k1')
    test_obj.set_working_path(self.path)
    test_obj.process_pdb_and_cif_files('tttt')
    test_obj.refinement_loop(
    n_macro_cycle        = 5,
    r_work_target        = 0.0001,
    print_during_refinement = True,
    sites                = False,
    u_iso                = False,
    transformations      = True)
    print test_obj
    print '='*173

  def test_refine_witout_transformations(self):
    """
    Test transformation refinement
    """
    print 'Running ',sys._getframe().f_code.co_name
    test_obj = ncs_refine_test(
      use_geometry_restraints=False,  # change to True
      real_data=False,
      sf_and_grads_algorithm = 'direct',
      target_name = 'ls_wunit_k1')
    test_obj.set_working_path(self.path)
    test_obj.process_pdb_and_cif_files('nnnn')
    test_obj.refinement_loop(
    n_macro_cycle        = 5,
    r_work_target        = 0.0001,
    print_during_refinement = True,
    sites                = True,
    u_iso                = True,
    transformations      = False)
    print test_obj
    print '='*173

  def test_altrnating_method(self):
    """
    Verify that alternating refinement method works properly
    """
    print 'Running ',sys._getframe().f_code.co_name
    test_obj = ncs_refine_test(
      use_geometry_restraints=False,  # change to True
      real_data=False,
      sf_and_grads_algorithm = 'direct',
      target_name = 'ls_wunit_k1')
    # test all three
    method = test_obj.use_refinement_method(
      sites=True,u_iso=True, transformations=True)
    for i in range(10):
      n = i % 3
      method = test_obj.iterate_refine_method(method,i)
      temp = [ v[2] for (k,v) in method.iteritems()]
      assert temp.count(True) == 1
      if n == 0: assert method['sites'][2] == True
      if n == 1: assert method['u_iso'][2] == True
      if n == 2: assert method['transformations'][2] == True
    # test just one method
    method = test_obj.use_refinement_method(
      sites=False,u_iso=True, transformations=False)
    for i in range(10):
      method = test_obj.iterate_refine_method(method,i)
      temp = [ v[2] for (k,v) in method.iteritems()]
      assert temp.count(True) == 1
      assert method['sites'][2] == False
      assert method['u_iso'][2] == True
      assert method['transformations'][2] == False
    # test two
    # test just one method
    method = test_obj.use_refinement_method(
      sites=False,u_iso=True, transformations=True)
    for i in range(10):
      n = i % 2
      method = test_obj.iterate_refine_method(method,i)
      temp = [ v[2] for (k,v) in method.iteritems()]
      assert temp.count(True) == 1
      assert method['sites'][2] == False
      if n == 0: assert method['u_iso'][2] == True
      if n == 1: assert method['transformations'][2] == True



  def tearDown(self):
    '''remove temp files and folder'''
    print 'Running ',sys._getframe().f_code.co_name
    os.chdir(self.currnet_dir)

def run():
  test_case = test_ncs_refinement()
  test_case.setUp()
  test_case.test_sites()
  test_case.test_adp()
  test_case.test_alternate()
  test_case.test_with_geometry_restraints()
  test_case.test_transformation_refinement()
  test_case.test_refine_witout_transformations()
  test_case.test_altrnating_method()
  test_case.tearDown()

if __name__=='__main__':
  run()
#
#   # Analyze code
#   import cProfile
#   import pstats
#   cProfile.run("run()",filename='cProfile_log')
#   p = pstats.Stats('cProfile_log')
#   p.sort_stats('time').print_stats(15)

