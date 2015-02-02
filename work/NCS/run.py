from __future__ import division
from iotbx.pdb.multimer_reconstruction import multimer
import mmtbx.refinement.minimization_ncs_constraints
from libtbx.test_utils import approx_equal
from scitbx.array_family import flex
from cctbx import xray
import mmtbx.f_model
import mmtbx.utils
import iotbx.pdb
import sys
from libtbx.utils import null_out
from iotbx import reflection_file_utils
import mmtbx.restraints

def process_pdb_file(pdb_file_name):
  pdb_inp_one_ncs = iotbx.pdb.input(file_name=pdb_file_name)
  cs = pdb_inp_one_ncs.crystal_symmetry_from_cryst1()
  m = multimer(
    file_name=pdb_file_name,
    round_coordinates=False,
    reconstruction_type='cau',
    error_handle=True,
    eps=1e-2)
  pdb_hierarchy_asu = m.assembled_multimer
  processed_pdb_files_srv = mmtbx.utils.process_pdb_file_srv(
    crystal_symmetry      = cs,
    stop_for_unknowns     = True,
    log                   = sys.stdout,
    use_neutron_distances = False)
  processed_pdb_file, junk = processed_pdb_files_srv.\
    process_pdb_files(
      raw_records = pdb_hierarchy_asu.as_pdb_string().splitlines())
  assert len(m.rotation_matrices)>0
  ncs_selection = flex.bool(
    pdb_inp_one_ncs.xray_structure_simple().scatterers().size(), True)
  return processed_pdb_file, m, ncs_selection

def get_fmodel(xray_structure, hkl_file):
  rfs = reflection_file_utils.reflection_file_server(
    crystal_symmetry = xray_structure.crystal_symmetry(),
    force_symmetry   = True,
    reflection_files = hkl_file,
    err              = null_out())
  determine_data_and_flags_result = mmtbx.utils.determine_data_and_flags(
    reflection_file_server  = rfs,
    keep_going              = True,
    log                     = null_out())
  f_obs = determine_data_and_flags_result.f_obs
  assert xray_structure.crystal_symmetry().is_similar_symmetry(
    f_obs.crystal_symmetry())
  r_free_flags = determine_data_and_flags_result.r_free_flags
  assert r_free_flags is not None
  print "f-obs labels:", f_obs.info().labels
  print "r-free-flags labels", r_free_flags.info().labels
  fmodel = mmtbx.f_model.manager(
    xray_structure = xray_structure,
    f_obs          = f_obs,
    r_free_flags   = r_free_flags,
    target_name    = "ml")
  fmodel.update_all_scales(update_f_part1_for=False)
  print "r_work, r_free: %6.4f %6.4f" % (fmodel.r_work(), fmodel.r_free())
  return fmodel

def get_weight(fmodel, grm):
  fmdc = fmodel.deep_copy()
  fmdc.xray_structure.shake_sites_in_place(mean_distance=0.3)
  fmdc.update_xray_structure(xray_structure = fmdc.xray_structure,
    update_f_calc=True)
  fmdc.xray_structure.scatterers().flags_set_grads(state=False)
  xray.set_scatterer_grad_flags(
    scatterers = fmdc.xray_structure.scatterers(),
    site       = True)
  # fmodel gradients
  gxc = flex.vec3_double(fmdc.one_time_gradients_wrt_atomic_parameters(
    site = True).packed())
  # manager restraints, energy sites gradients
  gc = grm.energies_sites(
    sites_cart        = fmdc.xray_structure.sites_cart(),
    compute_gradients = True).gradients
  gc_norm  = gc.norm()
  gxc_norm = gxc.norm()
  weight = 1.
  if(gxc_norm != 0.0):
    weight = gc_norm / gxc_norm
  return weight

def exercise(args):
  processed_args = mmtbx.utils.process_command_line_args(
    args=args, log=null_out())
  ppf, multimer_obj, ncs_selection = process_pdb_file(
    pdb_file_name = processed_args.pdb_file_names[0])
  pdb_hierarchy = ppf.all_chain_proxies.pdb_hierarchy
  xray_structure = ppf.xray_structure()
  ncs_selection = flex.bool(xray_structure.scatterers().size(),
    ncs_selection.iselection())
  fmodel = get_fmodel(
    xray_structure = xray_structure,
    hkl_file=processed_args.reflection_files)
  geometry = ppf.geometry_restraints_manager(
    show_energies                = False,
    plain_pairs_radius           = 5,
    assume_hydrogens_all_missing = True)
  restraints_manager = mmtbx.restraints.manager(
    geometry      = geometry,
    normalization = True)
  restraints_manager.crystal_symmetry = xray_structure.crystal_symmetry()
  #
  print "start r_factor: %6.4f" % fmodel.r_work()
  for macro_cycle in xrange(10):
    data_weight = get_weight(fmodel=fmodel, grm=restraints_manager)
    minimized = mmtbx.refinement.minimization_ncs_constraints.lbfgs(
      fmodel                       = fmodel,
      rotation_matrices            = multimer_obj.rotation_matrices,
      translation_vectors          = multimer_obj.translation_vectors,
      ncs_atom_selection           = ncs_selection,
      finite_grad_differences_test = False,
      geometry_restraints_manager  = restraints_manager,
      data_weight                  = data_weight,
      refine_sites                 = True,
      refine_u_iso                 = False)
    print fmodel.r_work(), data_weight
    assert approx_equal(fmodel.r_work(), minimized.fmodel.r_work())
    fmodel.update_all_scales(update_f_part1_for=False)

if(__name__ == "__main__"):
  exercise(sys.argv[1:])
