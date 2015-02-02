from cctbx.array_family import flex
import mmtbx.model
from mmtbx import monomer_library
import mmtbx.monomer_library.server
import mmtbx.monomer_library.pdb_interpretation
from mmtbx import utils
from cctbx import xray

pdb_str="""\n
CRYST1   10.000   10.000   10.000  90.00  90.00  90.00 P 21 21 21    4
ATOM      1  N   LYS     1       5.000   5.000   5.000  1.00 20.00           N
ATOM      1  N   LYS     2       6.000   5.000   5.000  1.00 20.00           N
ATOM      1  N   LYS     4       5.000   5.500   5.500  1.00 20.00           N
TER
END
"""

pdb_str2='''\n
CRYST1   20.000   20.000   20.000  90.00  90.00  90.00 P 1
SCALE1      0.050000  0.000000  0.000000        0.00000
SCALE2      0.000000  0.050000  0.000000        0.00000
SCALE3      0.000000  0.000000  0.050000        0.00000
ATOM      1  N   LYS     1       5.000   5.000   5.000  1.00 20.00           N
ATOM      2  N   LYS     2       6.000   5.000   5.000  1.00 20.00           N
ATOM      3  N   LYS     4       5.000   5.500   5.500  1.00 20.00           N
TER
HETATM    4  O   HOH     1       5.500   5.000   5.000  1.00 10.00           O
HETATM    5  O   HOH     2       6.500   5.000   5.000  1.00 10.00           O
HETATM    6  O   HOH     3       5.500   5.500   5.500  1.00 10.00           O
TER
END'''

def exercise():
  mon_lib_srv = monomer_library.server.server()
  ener_lib = monomer_library.server.ener_lib()
  processed_pdb_file = monomer_library.pdb_interpretation.process(
    mon_lib_srv    = mon_lib_srv,
    ener_lib       = ener_lib,
    raw_records    = pdb_str,
    force_symmetry = True)
  geometry = processed_pdb_file.geometry_restraints_manager(
    show_energies      = False,
    plain_pairs_radius = 5.0)
  xrs = processed_pdb_file.xray_structure()
  obj = geometry.get_nonbonded_clashscore(
    sites_cart  = xrs.sites_cart(),
    site_labels = xrs.scatterers().extract_labels(),
    hd_sel      = xrs.hd_selection())
  print obj.nb_clashscore_all_clashes
  print obj.nb_clashscore_due_to_sym_op
  print obj.nb_clashscore_without_sym_op
  for rec in obj.nb_clash_proxies_all_clashes:
    print 'out: ',rec
  print '*'*60
  restraints_manager = mmtbx.restraints.manager(geometry = geometry,
    normalization = False)
  mol = mmtbx.model.manager(
    restraints_manager = restraints_manager,
    xray_structure = xrs,
    pdb_hierarchy = processed_pdb_file.all_chain_proxies.pdb_hierarchy)
  #
  new_scatterers = flex.xray_scatterer(
    xrs.scatterers().size(),
    xray.scatterer(occupancy = 1, b = 10, scattering_type = "O"))
  new_sites_frac = xrs.unit_cell().fractionalize(xrs.sites_cart()+[0.5,0,0])
  new_scatterers.set_sites(new_sites_frac)
  new_xrs = xray.structure(
    special_position_settings = xrs,
    scatterers                = new_scatterers)
  mol.add_solvent(
    solvent_xray_structure = new_xrs,
    refine_occupancies     = False,
    refine_adp             = "isotropic")
  obj = mol.restraints_manager.geometry.get_nonbonded_clashscore(
    sites_cart  = mol.xray_structure.sites_cart(),
    site_labels = mol.xray_structure.scatterers().extract_labels(),
    hd_sel      = mol.xray_structure.hd_selection())
  print obj.nb_clashscore_all_clashes
  print obj.nb_clashscore_due_to_sym_op
  print obj.nb_clashscore_without_sym_op
  for rec in obj.nb_clash_proxies_all_clashes:
    print rec
  print '*'*60
  of = open("junk.pdb","w")
  mol.write_pdb_file(out=of)
  of.close()
  
  mon_lib_srv = monomer_library.server.server()
  ener_lib = monomer_library.server.ener_lib()
  processed_pdb_file = monomer_library.pdb_interpretation.process(
    mon_lib_srv    = mon_lib_srv,
    ener_lib       = ener_lib,
    raw_records    = pdb_str2,
    force_symmetry = True)
  geometry = processed_pdb_file.geometry_restraints_manager(
    show_energies      = False,
    plain_pairs_radius = 5.0)
  xrs = processed_pdb_file.xray_structure()
  obj = geometry.get_nonbonded_clashscore(
    sites_cart  = xrs.sites_cart(),
    site_labels = xrs.scatterers().extract_labels(),
    hd_sel      = xrs.hd_selection())
  print obj.nb_clashscore_all_clashes
  print obj.nb_clashscore_due_to_sym_op
  print obj.nb_clashscore_without_sym_op
  for rec in obj.nb_clash_proxies_all_clashes:
    print 'out: ',rec
  


if (__name__ == "__main__"):
  exercise()
