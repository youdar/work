from __future__ import division
import mmtbx.monomer_library.pdb_interpretation
import os

os.chdir(r'C:\Phenix\Dev\Work\work\NCS\junk\pdb_test')

pdb_processed_file = mmtbx.monomer_library.pdb_interpretation.run(args=[
  'crystall_asymmetric_unit_1tnv.pdb'],
  assume_hydrogens_all_missing=False,
  hard_minimum_nonbonded_distance=0.0,
  nonbonded_distance_threshold=None,
  substitute_non_crystallographic_unit_cell_if_necessary=True)

grm = pdb_processed_file.geometry_restraints_manager()

xrs = pdb_processed_file.xray_structure()
sites_cart = xrs.sites_cart()
site_labels = xrs.scatterers().extract_labels()
hd_sel = xrs.hd_selection()

bonded_pairs = grm.pair_proxies().bond_proxies

nb_clash_info = grm.get_nonbonded_clashscore(
  sites_cart=sites_cart,
  site_labels=site_labels,
  hd_sel=hd_sel)
