from libtbx import easy_pickle
from mmtbx import monomer_library
import mmtbx.monomer_library.server
import mmtbx.monomer_library.pdb_interpretation
from scitbx.array_family import flex

pdb_str0="""\n
CRYST1   10.000   10.000   10.000  90.00  90.00  90.00 P 21 21 21    4
ATOM      1  N   LYS     1       5.000   5.000   5.000  1.00 20.00           N
ATOM      1  N   LYS     2       6.000   5.000   5.000  1.00 20.00           N
ATOM      1  N   LYS     4       5.000   5.500   5.500  1.00 20.00           N
TER
END
"""

pdb_str1="""\n
CRYST1   10.000   10.000   10.000  90.00  90.00  90.00 P 21 21 21    4
ATOM      1  N   LYS     1       5.000   5.000   5.000  1.00 20.00           N
ATOM      1  N   LYS     2       5.100   5.000   5.000  1.00 20.00           N
ATOM      1  N   LYS     4       5.000   5.500   5.500  1.00 20.00           N
TER
END
"""

def exercise(pdb_str, i):
  print "-"*79
  of = open("m%s.pdb"%str(i), "w")
  print >> of, pdb_str
  of.close()
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
  assert xrs is not None
  # tst 1
  obj = geometry.get_nonbonded_clashscore(
    sites_cart  = xrs.sites_cart(),
    site_labels = xrs.scatterers().extract_labels(),
    hd_sel      = xrs.hd_selection())
  print obj.nb_clashscore_all_clashes
  print obj.nb_clashscore_due_to_sym_op
  print obj.nb_clashscore_without_sym_op
  print
  # tst 2
  sel = flex.bool([True, True, False])
  xrs = xrs.select(sel)
  geometry = geometry.select(selection=sel)
  obj = geometry.get_nonbonded_clashscore(
    sites_cart  = xrs.sites_cart(),
    site_labels = xrs.scatterers().extract_labels(),
    hd_sel      = xrs.hd_selection())
  print obj.nb_clashscore_all_clashes
  print obj.nb_clashscore_due_to_sym_op
  print obj.nb_clashscore_without_sym_op
  
if (__name__ == "__main__"):
  for i, pdb_str in enumerate([pdb_str0, pdb_str1]):
    exercise(pdb_str=pdb_str, i=i)

