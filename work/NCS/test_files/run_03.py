from __future__ import division
from iotbx.pdb.multimer_reconstruction import multimer
import iotbx.pdb
import mmtbx.f_model
from cctbx import xray
import scitbx.lbfgs
import getpass
import os
import sys
from scitbx.array_family import flex

ncs_1_copy="""\
MTRIX1   1  1.000000  0.000000  0.000000        0.00000    1
MTRIX2   1  0.000000  1.000000  0.000000        0.00000    1
MTRIX3   1  0.000000  0.000000  1.000000        0.00000    1
MTRIX1   2  0.496590 -0.643597  0.582393        0.00000
MTRIX2   2  0.867925  0.376088 -0.324443        0.00000
MTRIX3   2 -0.010221  0.666588  0.745356        0.00000
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

class minimizer(object):
  def __init__(self,
        fmodel,
        ncs_transformations_object=None,
        ncs_atom_selection = None,
        run_finite_grad_differences_test = False,
        max_iterations=100,
        sites = False,
        u_iso = False):
    """Implementing strict NCS to refinement minimization

    Arguments:
    fmodel : fmodel of the complete ASU
    ncs_transformation_object : information on the NCS to ASU
       transformations and chains. A multimer object
    ncs_atom_selection : boolean array for selection of atoms in the NCS.
       A flex bool array
    """
    self.fmodel = fmodel
    self.fmodel.xray_structure.scatterers().flags_set_grads(state=False)
    self.x_target_functor = self.fmodel.target_functor()
    self.sites = sites
    self.u_iso = u_iso
    self.ncs_to_asu = ncs_transformations_object
    self.run_finite_grad_differences_test = run_finite_grad_differences_test
    if run_finite_grad_differences_test:
      # perform gradient calc test
      self.buffer_max_grad = flex.double()
      self.buffer_calc_grad = flex.double()
    # xray structure of NCS chains for self.x
    ncs_fmodel_xrs = self.fmodel.xray_structure.select(ncs_atom_selection)
    if(self.sites):
      self.x = ncs_fmodel_xrs.sites_cart().as_double()
    if(self.u_iso):
      assert ncs_fmodel_xrs.scatterers().size() == \
        ncs_fmodel_xrs.use_u_iso().count(True)
      self.x = ncs_fmodel_xrs.extract_u_iso_or_u_equiv()
    # Use all scatterers for gradient calculations
    if(self.sites):
      xray.set_scatterer_grad_flags(
        scatterers = self.fmodel.xray_structure.scatterers(),
        site       = True)
    if(self.u_iso):
      sel = flex.bool(
        self.fmodel.xray_structure.scatterers().size(), True).iselection()
      self.fmodel.xray_structure.scatterers().flags_set_grad_u_iso(
        iselection = sel)
    self.minimizer = scitbx.lbfgs.run(
      target_evaluator=self,
      termination_params=scitbx.lbfgs.termination_parameters(
        max_iterations=max_iterations),
      exception_handling_params=scitbx.lbfgs.exception_handling_parameters(
        ignore_line_search_failed_rounding_errors=True,
        ignore_line_search_failed_step_at_lower_bound=True,
        ignore_line_search_failed_maxfev=True))
    self.fmodel.xray_structure.tidy_us()
    self.fmodel.xray_structure.apply_symmetry_sites()
    self.fmodel.update_xray_structure(
      xray_structure = self.fmodel.xray_structure,
      update_f_calc  = True)
    self.tested = 0
    if run_finite_grad_differences_test:
      if self.buffer_max_grad:
        print 'compare max_grad to calc_grad'
        for a,f in zip(self.buffer_max_grad, self.buffer_calc_grad):
          print '{0:10.5f}   {1:10.5f}  delta = {2:10.5f}'.format(a,f,abs(a-f))
        print '-'*45
        diff = flex.abs(self.buffer_max_grad - self.buffer_calc_grad)
        s = diff < 1.e-3
        if(s.size()>0 and s.count(True)*100./s.size()>50):
          self.tested += 1

  def compute_functional_and_gradients(self,compute_gradients=True):
    """(bool) -> float, flex.double array
    Function which calculates the target function and gradients.
    It is called by the lbfgs minimizer

    Argument:
    compute_gradients : When True gradients are calculated
    """
    if(self.sites):
      self.update_model_sites()
    elif(self.u_iso):
      self.update_model_asu_b_factors()
    self.fmodel.update_xray_structure(
      xray_structure = self.fmodel.xray_structure,
      update_f_calc  = True)
    tgx = self.x_target_functor(compute_gradients=compute_gradients)
    if(self.sites):
      tx = tgx.target_work()
      f = tx
      if compute_gradients:
        gx = flex.vec3_double(tgx.\
          gradients_wrt_atomic_parameters(site=True).packed())
        g = self.average_grad(grad=gx,apply_rotation=True).as_double()
        if self.run_finite_grad_differences_test:
          self.finite_difference_test(g)
    if(self.u_iso):
      tx = tgx.target_work()
      f = tx
      if compute_gradients:
        gx = tgx.gradients_wrt_atomic_parameters(u_iso=True)
        g = self.average_grad(grad=gx,apply_rotation=False).as_double()
    if not compute_gradients:
      g = None
    return f, g

  def update_model_sites(self,x=None):
    """
    update fmodel using a complete ASU

    Argument
    x : sites coordinates of a single NCS
    """
    if not x:
      x = self.x
    # rebuild the complete ASU coordinate from the NCS
    x_asu = self.rebuild_asu_from_ncs_coordinates(x)
    self.fmodel.xray_structure.set_sites_cart(
        sites_cart = flex.vec3_double(x_asu))

  def update_model_asu_b_factors(self,x_asu=None):
    """
    update fmodel using a complete set of B-factors
    All B-factors are the same

    Argument
    x_asu : B-factors of a single NCS
    """
    if not x_asu:
      x_asu = self.rebuild_asu_b_factors()
    self.fmodel.xray_structure.set_u_iso(values = x_asu)

  def rebuild_asu_b_factors(self):
    """
    """
    n = self.ncs_to_asu.number_of_transforms + 1
    x = list(self.x) * n
    return flex.double(x)

  def rebuild_asu_from_ncs_coordinates(self, x):
    """ apply rotation and translation to x

    Argument:
    x : sites coordinates of single NCS

    returns:
    new_x : coordinates of the original x and all coordinates resulting
       from application of rotation and translation.
       type scitbx_array_family_flex_ext.double
    """
    rotations = self.ncs_to_asu.rotation_matrices
    translations =  self.ncs_to_asu.translation_vectors
    assert len(rotations)==len(translations)
    new_x = list(x)
    x = flex.vec3_double(x)
    for r,t in zip(rotations,translations):
      tmp_x = r.elems*x + t
      new_x += list(tmp_x.as_double())
    return flex.double(new_x)

  def average_grad(self,grad,apply_rotation=False):
    """(vec3_double,bool) -> vec3_double

    Argument:
    grad : the gradient of the complete ASU
    apply_rotation : If true, apply NCS rotation before averaging

    Returns:
    g_ave : The average the gradients of all NCS copies in the ASU
    """
    n = self.ncs_to_asu.number_of_transforms
    # gradients of the first NCS copy
    ncs_end = len(grad)//(n+1)
    assert ncs_end*(n+1)==len(grad)
    g_ave = grad[:ncs_end]
    for i in range(n):
      g = grad[ncs_end*(i+1):ncs_end*(i+2)]
      if apply_rotation:
        # multiply the transpose of the rotation of each NCS copy
        # gradients, by the gradients
        rt = self.ncs_to_asu.rotation_matrices[i].transpose().elems
        g = rt*g
      g_ave += g
    # average the NCS copies contributions
    g_ave = g_ave.as_double()/(n+1)
    if apply_rotation: g_ave = flex.vec3_double(g_ave)
    assert type(grad)==type(g_ave)
    return g_ave

  def finite_difference_test(self,g):
    """
    Run basic gradient test. compare numerical estimate gradient to
    the largest calculated one. using t'(x)=(t(x+d)-t(x-d))/(2d)

    Argument:
     g : gradient, flex array
    """
    if(self.fmodel.r_work()>1.e-3):
      g = g.as_double()
      d = 1.e-5
      # find the index of the max gradient value
      i_g_max = flex.max_index(flex.abs(g))
      x_d = self.x
      # calc t(x+d)
      x_d[i_g_max] = self.x[i_g_max] + d
      self.update_model_sites(x = x_d)
      self.fmodel.update_xray_structure(update_f_calc=True)
      t1,_ = self.compute_functional_and_gradients(compute_gradients=False)
      # calc t(x-d)
      x_d[i_g_max] = self.x[i_g_max] - d
      self.update_model_sites(x = x_d)
      del x_d
      self.fmodel.update_xray_structure(update_f_calc=True)
      t2,_ = self.compute_functional_and_gradients(compute_gradients=False)
      # Return fmodel to the correct coordinates values
      self.update_model_sites(x = self.x)
      self.fmodel.update_xray_structure(update_f_calc=True)
      self.buffer_max_grad.append(g[i_g_max])
      self.buffer_calc_grad.append((t1-t2)/(d*2))

def save_pdb_file(macro_cycle,fmodel,m_shaken,u_iso,sites):
  """
  Save pdb file for visualization
  """
  method = ''
  if sites: method += '_sites'
  if u_iso: method += '_u_iso'
  # fn = 'refinement_by{0}_step_{1}.pdb'.format(method,macro_cycle)
  n = str(macro_cycle)
  if len(n)==1: n = '0'+n
  fn = 'refinement{}.pdb'.format(n)
  xrs_refined = fmodel.xray_structure
  m_shaken.assembled_multimer.adopt_xray_structure(xrs_refined)
  m_shaken.write(fn)

def create_pymol_movie():
  """create pymol movie
  """
  from glob import glob
  import pymol
  file_list = glob("refinement*.pdb")
  pymol.finish_launching()
  pymol.cmd.bg_color('white')
  for i,fn in enumerate(file_list):
    # pymol.cmd.load(fn,"mov",state=i)
    # pymol.cmd.load('full_asu.pdb',"mov",state=i)
    pymol.cmd.load(fn)
    pymol.cmd.load('full_asu.pdb')
    pymol.cmd.frame(1)
    pymol.cmd.mview('store')
  #
  pymol.cmd.mset("1 -%d"%len(file_list))

def run(
        n_macro_cycle=10,
        sites=True,
        u_iso=False,
        run_finite_grad_differences_test = False):
  """
  Arguments:
  __________
  n_macro_cycle : Number of refinement macro cycles
  """
  # 1 NCS copy: starting template to generate whole asu; place into P1 box
  pdb_inp = iotbx.pdb.input(source_info=None, lines=ncs_1_copy)
  mtrix_object = pdb_inp.process_mtrix_records()
  ph = pdb_inp.construct_hierarchy()
  xrs = pdb_inp.xray_structure_simple()
  xrs_one_ncs = xrs.orthorhombic_unit_cell_around_centered_scatterers(
    buffer_size=8)
  ph.adopt_xray_structure(xrs_one_ncs)
  of = open("one_ncs_in_asu.pdb", "w")
  print >> of, mtrix_object.format_MTRIX_pdb_string()
  print >> of, ph.as_pdb_string(crystal_symmetry=xrs_one_ncs.crystal_symmetry())
  of.close()
  # 1 NCS copy -> full asu (expand NCS). This is the answer-strucure
  m = multimer("one_ncs_in_asu.pdb",'cau',error_handle=True,eps=1e-2)
  assert m.number_of_transforms == 2, m.number_of_transforms
  xrs_asu = m.assembled_multimer.extract_xray_structure(
    crystal_symmetry = xrs_one_ncs.crystal_symmetry())
  m.write("full_asu.pdb")
  assert xrs_asu.crystal_symmetry().is_similar_symmetry(
    xrs_one_ncs.crystal_symmetry())
  # Generate Fobs from answer structure
  f_obs = abs(xrs_asu.structure_factors(d_min=2, algorithm="direct").f_calc())
  r_free_flags = f_obs.generate_r_free_flags()
  mtz_dataset = f_obs.as_mtz_dataset(column_root_label="F-obs")
  mtz_dataset.add_miller_array(
    miller_array=r_free_flags,
    column_root_label="R-free-flags")
  mtz_object = mtz_dataset.mtz_object()
  mtz_object.write(file_name = "data.mtz")
  # Shake structure - subject to refinement input
  xrs_shaken = xrs_one_ncs.deep_copy_scatterers()
  if sites:
    xrs_shaken.shake_sites_in_place(mean_distance=0.3)
  if u_iso:
    xrs_shaken.shake_adp()
  ph.adopt_xray_structure(xrs_shaken)
  of = open("one_ncs_in_asu_shaken.pdb", "w")
  print >> of, mtrix_object.format_MTRIX_pdb_string()
  print >> of, ph.as_pdb_string(crystal_symmetry=xrs.crystal_symmetry())
  of.close()
  ### Refinement
  params = mmtbx.f_model.sf_and_grads_accuracy_master_params.extract()
  params.algorithm = "direct"
  # Get the xray_structure of the shaken ASU
  m_shaken = multimer(
    pdb_input_file_name="one_ncs_in_asu_shaken.pdb",
    reconstruction_type='cau',error_handle=True,eps=1e-2)
  xrs_shaken_asu = m_shaken.assembled_multimer.as_pdb_input().\
    xray_structure_simple(crystal_symmetry=xrs_one_ncs.crystal_symmetry())
  # Save the shaken ASU for inspection
  m_shaken.write(pdb_output_file_name='asu_shaken.pdb')
  # Create a boolean selection string for selecting chains in NCS
  selection_str = 'chain A'
  ncs_selection = m_shaken.assembled_multimer.\
    atom_selection_cache().selection(selection_str)
  fmodel = mmtbx.f_model.manager(
    f_obs                        = f_obs,
    r_free_flags                 = r_free_flags,
    xray_structure               = xrs_shaken_asu,
    sf_and_grads_accuracy_params = params,
    target_name                  = "ls_wunit_k1")
  print "start r_factor: %6.4f" % fmodel.r_work()
  refine_method = 'sites'
  for macro_cycle in xrange(n_macro_cycle):
    # refine coordinates
    if(sites):
      minimized = minimizer(
        fmodel = fmodel,
        ncs_transformations_object=m,
        ncs_atom_selection = ncs_selection,
        run_finite_grad_differences_test = run_finite_grad_differences_test,
        sites = True)
      print "  macro_cycle %3d (sites) r_factor: %6.4f"%(macro_cycle,
        fmodel.r_work())
    # refine ADPs
    if(u_iso):
      minimized = minimizer(
        fmodel = fmodel,
        ncs_transformations_object=m,
        ncs_atom_selection = ncs_selection,
        run_finite_grad_differences_test = run_finite_grad_differences_test,
        u_iso = True)
      print "  macro_cycle %3d (adp)   r_factor: %6.4f"%(macro_cycle,
        fmodel.r_work())
    if (0): save_pdb_file(
      macro_cycle=macro_cycle,
      fmodel=fmodel,
      m_shaken=m_shaken,
      u_iso=u_iso,
      sites=sites)
  if (1): save_pdb_file(
    macro_cycle=macro_cycle,
    fmodel=fmodel,
    m_shaken=m_shaken,
    u_iso=u_iso,
    sites=sites)
  # create_pymol_movie()

def set_test_folder():
    """
    Change working directory to avoid littering of
    phenix_sources\phenix_regression\development\ncs_constraints.py
    """
    username = getpass.getuser()
    if username.lower() == 'youval':
      osType = sys.platform
      if osType.startswith('win'):
        tempdir = (r'C:\Phenix\Dev\Work\work\NCS\junk')
      else:
        tempdir = ('/net/cci/youval/Work/work/NCS/junk')
      os.chdir(tempdir)

if __name__ == "__main__":
  set_test_folder()
  run(n_macro_cycle=40,
      sites=True,
      u_iso=False,
      run_finite_grad_differences_test = False)
