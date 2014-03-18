from __future__ import division
from iotbx.pdb.multimer_reconstruction import multimer
from scitbx.array_family import flex
from libtbx import group_args
import mmtbx.f_model
import scitbx.lbfgs
from iotbx import pdb
from cctbx import xray
import sys,os

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
      # Coordiantes
      self.x = self.fmodel.xray_structure.sites_cart().as_double()
    if(self.u_iso):
      # anisotropic displacement factors
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
    """() -> float,flex.double array

    This methode is called from:
    phenix_sources\cctbx_project\scitbx\lbfgs\__init__.py

    It returns the target_work and gradients to lbfgs
    """
    print '='*80
    if self.sites:
      self.fmodel.xray_structure.set_sites_cart(
        sites_cart = flex.vec3_double(self.x))
      print len(flex.vec3_double(self.x))
    if self.u_iso:
      self.fmodel.xray_structure.set_u_iso(values = self.x)
    self.fmodel.update_xray_structure(
      xray_structure = self.fmodel.xray_structure,
      update_f_calc  = True)
    tgx = self.x_target_functor(compute_gradients=True)
    if self.sites:
      tx = tgx.target_work()
      gx = flex.vec3_double(
        tgx.gradients_wrt_atomic_parameters(site=True).packed())
      f = tx
      g = gx
    if self.u_iso:
      tx = tgx.target_work()
      gx = tgx.gradients_wrt_atomic_parameters(u_iso=True)
      f = tx
      g = gx
    # When we have MTRIX records, use only the
    # gradients of the first NCS copy
    # ncs_end = len(g)//(self.n_ncs_mtrix+1)
    # assert ncs_end*(self.n_ncs_mtrix+1)==len(g)
    # g = g[:ncs_end]
    return f, g.as_double()

def get_inputs(file_name):
  '''
  obj = get_inputs(pdb_file_name)

  Retruns:
  --------
  obj.pdb_hierarchy
  obj.xray_structure
  '''
  pdb_inp = pdb.input(file_name=file_name)
  return group_args(
    pdb_hierarchy  = pdb_inp.construct_hierarchy(),
    xray_structure = pdb_inp.xray_structure_simple())

def run(file_to_refine,f_obs,r_free_flags,n_macro_cycle=10,r_work_target=0.0):
  """
  Arguments:
  ---------
  file_to_refine: a shaken copy of one NCS copy
                  used to produce the reference ASU
  f_obs: observed frequency
  r_free_flags: R-free-flags, flags for cross-validation data
  n_macro_cycle: Number of refinement cycles
  r_work_target: Refinement will stop when r_work <= r_work_target
  """
  crystal_symmetry = f_obs.crystal_symmetry()
  # Data to refine
  m = multimer(
    pdb_input_file_name=file_to_refine,
    reconstruction_type='cau',error_handle=True,eps=1e-2)
  # Keep original crystal symetry
  xrs_poor = m.assembled_multimer.as_pdb_input().\
    xray_structure_simple(crystal_symmetry=crystal_symmetry)
  # get fmodel
  params = mmtbx.f_model.sf_and_grads_accuracy_master_params.extract()
  params.algorithm = "direct"
  fmodel = mmtbx.f_model.manager(
    f_obs                        = f_obs,
    r_free_flags                 = r_free_flags,
    xray_structure               = xrs_poor,
    sf_and_grads_accuracy_params = params,
    target_name                  = "ls_wunit_kunit")
  # refinement loop
  print "start r_factor: %6.4f" % fmodel.r_work()
  for macro_cycle in xrange(n_macro_cycle):
    if(1):
      # refine coordinates
      minimized = minimizer(
        fmodel = fmodel,
        n_ncs_mtrix = m.number_of_transforms,
        sites = True)
      print "  macro_cycle %3d (sites) r_factor: %6.4f"%(macro_cycle, fmodel.r_work())
    if(0):
      # refine ADPs (atomic displacement parameter)
      minimized = minimizer(
        fmodel = fmodel,
        n_ncs_mtrix = m.number_of_transforms,
        u_iso = True)
      print "  macro_cycle %3d (adp)   r_factor: %6.4f"%(macro_cycle, fmodel.r_work())
      if fmodel.r_work() <= r_work_target: break
  if(1):
    m.assembled_multimer.adopt_xray_structure(fmodel.xray_structure)
    m.write(pdb_output_file_name="refined.pdb")

if __name__=='__main__':
  osType = sys.platform
  if osType.startswith('win'):
    tempdir = r'C:\Phenix\Dev\Work\work\NCS\junk'
  else:
    tempdir = '/net/cci/youval/Work/work/NCS/junk'
  os.chdir(tempdir)
  # Test files
  file_to_refine = 'ncs1_shaken.pdb'
  file_reference = 'asu0.pdb'
  # get xray_structure from reference PDB file
  inp = get_inputs(file_name=file_reference)
  xrs = inp.xray_structure
  # simulate Fobs
  f_obs = abs(xrs.structure_factors(d_min=1.0, algorithm="direct").f_calc())
  r_free_flags = f_obs.generate_r_free_flags()
  # Start refinement
  run(file_to_refine=file_to_refine,
      n_macro_cycle=2,
      f_obs=f_obs,
      r_free_flags=r_free_flags)

