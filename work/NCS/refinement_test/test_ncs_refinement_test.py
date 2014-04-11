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
    n_macro_cycle=100,
    r_work_target=0.00001,
    sites=True,
    u_iso=False,
    alternate_refinement=False)
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
    n_macro_cycle=50,
    r_work_target=0.00001,
    sites=False,
    u_iso=True,
    alternate_refinement=False)
    print test_obj
    print '='*173

  def test_adp_and_site(self):
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
    n_macro_cycle=50,
    r_work_target=0.00001,
    sites=False,
    u_iso=True,
    alternate_refinement=True)
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
    n_macro_cycle=50,
    r_work_target=0.0001,
    sites=True,
    u_iso=False,
    alternate_refinement=False)
    print test_obj
    print '='*173

  def tearDown(self):
    '''remove temp files and folder'''
    print 'Running ',sys._getframe().f_code.co_name
    os.chdir(self.currnet_dir)

def run():
  test_case = test_ncs_refinement()
  test_case.setUp()
  test_case.test_sites()
  test_case.test_adp()
  test_case.test_adp_and_site()
  test_case.test_with_geometry_restraints()
  test_case.tearDown()

if __name__=='__main__':
  run()

