


from __future__ import division
import mmtbx.monomer_library.pdb_interpretation
import iotbx.mtz
from cctbx.array_family import flex
import time
from mmtbx import monomer_library
import mmtbx.refinement.real_space.fit_residues
import iotbx.pdb
from libtbx import group_args



def run():

  
  # Get pdb file info
  pdb_inp = iotbx.pdb.input(file_name='')
  xrs_answer = pdb_inp.xray_structure_simple()
  f_calc = xrs_answer.structure_factors(d_min = d_min).f_calc()
  fft_map = f_calc.fft_map(resolution_factor=resolution_factor)
  fft_map.apply_sigma_scaling()
  target_map = fft_map.real_map_unpadded()
  mtz_dataset = f_calc.as_mtz_dataset(column_root_label = "FCmap")
  mtz_object = mtz_dataset.mtz_object()
  mtz_object.write(file_name = "answer_%s.mtz"%prefix)
  
  
  # poor
  mon_lib_srv = monomer_library.server.server()
  processed_pdb_file = monomer_library.pdb_interpretation.process(
    mon_lib_srv              = mon_lib_srv,
    ener_lib                 = monomer_library.server.ener_lib(),
    raw_records              = flex.std_string(pdb_str_poor.splitlines()),
    strict_conflict_handling = True,
    force_symmetry           = True,
    log                      = None)
  pdb_hierarchy_poor = processed_pdb_file.all_chain_proxies.pdb_hierarchy
  xrs_poor = processed_pdb_file.xray_structure()
  sites_cart_poor = xrs_poor.sites_cart()
  pdb_hierarchy_poor.write_pdb_file(file_name = "poor_%s.pdb"%prefix)
  dist = xrs_answer.mean_distance(other = xrs_poor)
  
  target_map_object = group_args(
    data             = target_map,
    miller_array     = f_calc,
    crystal_gridding = fft_map)
  
  
  sm = mmtbx.refinement.real_space.structure_monitor(
    pdb_hierarchy               = pdb_hierarchy_poor,
    xray_structure              = xrs_poor,
    target_map_object           = target_map_object,
    geometry_restraints_manager = grm.geometry)

def find_nonbonded():
  '''find nonbonded clashing pairs of atoms'''
  # Get bonds info using the atoms coordinates (site_cart) xray structure
  bond_proxies_simple = self.geometry_restraints_manager.pair_proxies(
    sites_cart = self.xray_structure.sites_cart()).bond_proxies.simple
  # Collect bonded pairs
  bonded_i_seqs = []
  for bp in bond_proxies_simple:
    bonded_i_seqs.append(bp.i_seqs)
  pair_asu_table = self.xray_structure.pair_asu_table(
    distance_cutoff=self.clash_threshold)
  pair_sym_table = pair_asu_table.extract_pair_sym_table()
  atom_pairs_i_seqs = pair_sym_table.simple_edge_list()
  nonbonded_pairs = list(set(atom_pairs_i_seqs).difference(set(bonded_i_seqs)))
  
def set_working_path():
  # locate the directory containing the log files
  osType = sys.platform
  if osType.startswith('win'):
      directory_path = 'c:\Phenix\Dev\Work\work\Clashes\junk'
  else:
      directory_path = '/net/cci-filer2/raid1/home/youval/Work/work/Clashes/junk'
  os.chdir(directory_path)

if __name__=='__main__':
  # set work folder
  set_working_path()
  run()