from __future__ import division
from iotbx.pdb.multimer_reconstruction import multimer
from iotbx import pdb
import iotbx.pdb
import mmtbx.f_model
from cctbx import xray
import scitbx.lbfgs
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
        max_iterations=100,
        n_ncs_mtrix=0,
        sites = False,
        u_iso = False):
    self.fmodel = fmodel
    self.fmodel.xray_structure.scatterers().flags_set_grads(state=False)
    self.x_target_functor = self.fmodel.target_functor()
    self.n_ncs_mtrix = n_ncs_mtrix
    self.sites = sites
    self.u_iso = u_iso
    if(self.sites):
      self.x = self.fmodel.xray_structure.sites_cart().as_double()
    if(self.u_iso):
      assert self.fmodel.xray_structure.scatterers().size() == \
        self.fmodel.xray_structure.use_u_iso().count(True)
      self.x = self.fmodel.xray_structure.extract_u_iso_or_u_equiv()
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

  def compute_functional_and_gradients(self):
    if(self.sites):
      self.fmodel.xray_structure.set_sites_cart(
        sites_cart = flex.vec3_double(self.x))
    if(self.u_iso):
      self.fmodel.xray_structure.set_u_iso(values = self.x)
    self.fmodel.update_xray_structure(
      xray_structure = self.fmodel.xray_structure,
      update_f_calc  = True)
    tgx = self.x_target_functor(compute_gradients=True)
    if(self.sites):
      tx = tgx.target_work()
      gx = flex.vec3_double(tgx.\
        gradients_wrt_atomic_parameters(site=True).packed())
      f = tx
      g = gx
    if(self.u_iso):
      tx = tgx.target_work()
      gx = tgx.gradients_wrt_atomic_parameters(u_iso=True)
      f = tx
      g = gx
    # When we have MTRIX records, use only the
    # gradients of the first NCS copy
    ncs_end = len(g)//(self.n_ncs_mtrix+1)
    assert ncs_end*(self.n_ncs_mtrix+1)==len(g)
    g[:ncs_end]
    return f, g.as_double()

def run():
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
  xrs_shaken.shake_sites_in_place(mean_distance=0.3)
  ph.adopt_xray_structure(xrs_shaken)
  of = open("one_ncs_in_asu_shaken.pdb", "w")
  print >> of, mtrix_object.format_MTRIX_pdb_string()
  print >> of, ph.as_pdb_string(crystal_symmetry=xrs.crystal_symmetry())
  of.close()
  # Refinement
  params = mmtbx.f_model.sf_and_grads_accuracy_master_params.extract()
  params.algorithm = "direct"

  m_shaken = multimer(
    pdb_input_file_name="one_ncs_in_asu_shaken.pdb",
    reconstruction_type='cau',error_handle=True,eps=1e-2)
  xrs_shaken_asu = m_shaken.assembled_multimer.as_pdb_input().\
    xray_structure_simple(crystal_symmetry=xrs_one_ncs.crystal_symmetry())

  fmodel = mmtbx.f_model.manager(
    f_obs                        = f_obs,
    r_free_flags                 = r_free_flags,
    xray_structure               = xrs_shaken_asu,
    sf_and_grads_accuracy_params = params,
    target_name                  = "ls_wunit_k1")
  print "start r_factor: %6.4f" % fmodel.r_work()
  for macro_cycle in xrange(10):
    # refine coordinates
    if(1):
      minimized = minimizer(
        fmodel = fmodel,
        n_ncs_mtrix = m.number_of_transforms,
        sites = True)
      print "  macro_cycle %3d (sites) r_factor: %6.4f"%(macro_cycle,
        fmodel.r_work())
    # refine ADPs
    if(0):
      minimized = minimizer(
        fmodel = fmodel,
        n_ncs_mtrix = m.number_of_transforms,
        u_iso = True)
      print "  macro_cycle %3d (adp)   r_factor: %6.4f"%(macro_cycle,
        fmodel.r_work())

if __name__ == "__main__":
  run()

