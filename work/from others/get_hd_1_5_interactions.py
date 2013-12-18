
from __future__ import division
import sys

def run (args) :
  from mmtbx.monomer_library import pdb_interpretation
  from cctbx import crystal
  pdb = pdb_interpretation.run(args=args,
    substitute_non_crystallographic_unit_cell_if_necessary=True)
  geo = pdb.geometry_restraints_manager()
  xrs = pdb.xray_structure()
  sites_cart = xrs.sites_frac()
  scatterers = xrs.scatterers()
  hd_sel = xrs.hd_selection()
  table_bonds = geo.shell_sym_tables[0]
  table_1_4 = geo.shell_sym_tables[2]
  assert len(table_1_4) == len(sites_cart)
  for i_seq, sc in enumerate(scatterers) :
    print sc.label
    if (hd_sel[i_seq]) :
      bonded = table_bonds[i_seq].keys()
      for j_seq in bonded :
        for k_seq in table_1_4[j_seq].keys() :
          print "  H/D 1_5 interaction:", scatterers[k_seq].label
    else :
      for j_seq in table_1_4[i_seq].keys() :
        for k_seq in table_bonds[j_seq].keys() :
          if (hd_sel[k_seq]) :
            print "  H/D 1_5 interaction:", scatterers[k_seq].label

if (__name__ == "__main__") :
  run(sys.argv[1:])
