from __future__ import division
import mmtbx.monomer_library.pdb_interpretation as pdb_inter
import cctbx.geometry_restraints.nonbonded_overlaps as cs
from libtbx.utils import null_out
import iotbx.pdb

test_pdb_str = """\
CRYST1   44.060   35.400   48.340  90.00  95.00  90.00 C 1 2 1       4
ATOM      1  N   GLY A   1      -6.724   4.519  10.133  1.00 16.77           N
ATOM      2  CA  GLY A   1      -7.194   4.166   8.745  1.00 16.57           C
ATOM      3  C   GLY A   1      -6.271   3.120   8.177  1.00 16.16           C
ATOM      4  O   GLY A   1      -5.516   2.473   8.927  1.00 16.78           O
ATOM      5  N   ASN A   2      -6.301   2.953   6.856  1.00 15.02           N
ATOM      6  CA  ASN A   2      -5.313   2.093   6.179  1.00 14.10           C
ATOM      7  C   ASN A   2      -3.913   2.586   6.388  1.00 13.13           C
ATOM      8  O   ASN A   2      -3.663   3.793   6.355  1.00 11.91           O
ATOM      9  CB  ASN A   2      -5.585   1.992   4.699  1.00 15.38           C
ATOM     10  CG  ASN A   2      -6.959   1.462   4.424  1.00 14.08           C
ATOM     11  OD1 ASN A   2      -7.284   0.331   4.822  1.00 17.46           O
ATOM     12  ND2 ASN A   2      -7.807   2.298   3.816  1.00 11.72           N
ATOM     13  N   ASN A   3      -3.004   1.632   6.572  1.00 12.26           N
ATOM     14  CA  ASN A   3      -1.613   1.936   6.869  1.00 11.74           C
ATOM     15  C   ASN A   3      -0.635   1.404   5.819  1.00 11.10           C
ATOM     16  O   ASN A   3      -0.628   0.202   5.513  1.00 10.42           O
ATOM     17  CB  ASN A   3      -1.246   1.358   8.256  1.00 12.15           C
ATOM     18  CG  ASN A   3       0.193   1.704   8.681  1.00 12.82           C
ATOM     19  OD1 ASN A   3       0.544   2.886   8.837  1.00 15.05           O
ATOM     20  ND2 ASN A   3       1.027   0.674   8.850  1.00 13.48           N
ATOM     21  N   GLN A   4       0.718   2.275   5.803  1.00 10.29           N
ATOM     22  CA  GLN A   4       1.399   2.012   4.510  1.00 10.53           C
ATOM     23  C   GLN A   4       2.797   2.684   4.444  1.00 10.24           C
ATOM     24  O   GLN A   4       2.943   3.917   4.484  1.00  8.86           O
ATOM     25  CB  GLN A   4       0.545   2.406   3.297  1.00  9.80           C
ATOM     26  CG  GLN A   4       1.072   1.800   1.978  1.00 10.25           C
ATOM     27  CD  GLN A   4       0.565   2.510   0.717  1.00 12.43           C
ATOM     28  OE1 GLN A   4       0.710   3.745   0.575  1.00 14.62           O
ATOM     29  NE2 GLN A   4      -0.007   1.723  -0.230  1.00  9.05           N
ATOM     30  N   GLN A   5       3.828   1.858   4.418  1.00 10.38           N
ATOM     31  CA  GLN A   5       4.819   1.038   5.107  1.00 11.39           C
ATOM     32  C   GLN A   5       5.215  -0.166   4.260  1.00 11.52           C
ATOM     33  O   GLN A   5       4.376  -0.926   3.752  1.00 12.05           O
ATOM     34  CB  GLN A   5       4.342   0.629   6.505  1.00 11.96           C
ATOM     35  CG  GLN A   5       4.135   1.841   7.417  1.00 10.81           C
ATOM     36  CD  GLN A   5       3.241   1.514   8.568  1.00 13.10           C
ATOM     37  OE1 GLN A   5       2.035   1.354   8.386  1.00 10.65           O
ATOM     38  NE2 GLN A   5       3.822   1.429   9.781  1.00 12.30           N
TER
HETATM   61  O   HOH A   8      -0.511   4.797  12.393  1.00 22.62           O
HETATM   62  O   HOH A   9      -0.513   4.516  12.150  1.00 19.71           O
TER
HETATM 2800  S   SO4 A 701      -3.889   1.786  10.440  1.00 55.67           S
HETATM 2801  O1  SO4 A 701      -3.645   1.548   9.055  1.00 57.05           O
HETATM 2802  O2  SO4 A 701      -4.464   3.089  10.608  1.00 55.53           O
HETATM 2803  O3  SO4 A 701      -4.744   0.755  10.958  1.00 56.44           O
HETATM 2804  O4  SO4 A 701      -2.664   1.753  11.146  1.00 56.08           O
ATOM     60  O   HOH C  37      -0.639  -0.486   5.076  1.00  0.00           O
END
"""

def run():
  """
  Looking of the effect of adding hydrogen atoms on symmetry clashes

  Note that adding H changes the symmetry related clashscore
  """
  fn = 'test_file_no_h.pdb'
  open(fn,'w').write(test_pdb_str)
  clash_score_info_no_h = process_clash_score(fn)

  fn_with_h = 'test_file_with_h.pdb'
  pdb_with_h, h_were_added = cs.check_and_add_hydrogen(
        file_name=fn,
        allow_multiple_models=False,
        log=null_out())

  open(fn_with_h,'w').write(pdb_with_h)
  clash_score_info_with_h = process_clash_score(fn_with_h)

  clash_score_info_no_h.show()
  print '='*50
  clash_score_info_with_h.show()

  # compare the ca coordinate before and after adding H
  pdb_inp_no_h = iotbx.pdb.input(source_info=None, lines=test_pdb_str)
  pdb_inp_with_h = iotbx.pdb.input(source_info=None, lines=pdb_with_h)
  ph_no_h = pdb_inp_no_h.construct_hierarchy()
  ph_with_h = pdb_inp_with_h.construct_hierarchy()
  cashe_no_h = ph_no_h.atom_selection_cache().selection
  cashe_with_h = ph_with_h.atom_selection_cache().selection
  select_str = 'name ca or name o or name n'
  sites_cart_no_h = ph_no_h.atoms().select(cashe_no_h(select_str))
  sites_cart_with_h = ph_with_h.atoms().select(cashe_with_h(select_str))
  xyz_no_h = sites_cart_no_h.extract_xyz()
  xyz_with_h = sites_cart_with_h.extract_xyz()
  assert sites_cart_no_h.size() == sites_cart_with_h.size()
  #
  x = (xyz_no_h - xyz_with_h).as_double()
  print x.min_max_mean().as_tuple()






  # check f_calc

  #
  xrs_no_h = pdb_inp_no_h.xray_structure_simple()
  xrs_with_h = pdb_inp_with_h.xray_structure_simple()
  #
  f_calc_no_h = xrs_no_h.structure_factors(d_min = 2).f_calc()
  f_calc_with_h = xrs_with_h.structure_factors(d_min = 2).f_calc()


def process_clash_score(file_name):
  pdb_processed_file = pdb_inter.run(
      args=[file_name],
      assume_hydrogens_all_missing=False,
      hard_minimum_nonbonded_distance=0.0,
      nonbonded_distance_threshold=None,
      substitute_non_crystallographic_unit_cell_if_necessary=True,
      log=null_out())

  grm = pdb_processed_file.geometry_restraints_manager()
  xrs = pdb_processed_file.xray_structure()
  sites_cart = xrs.sites_cart()
  site_labels = xrs.scatterers().extract_labels()
  hd_sel = xrs.hd_selection()
  macro_mol_sel = cs.get_macro_mol_sel(pdb_processed_file)

  clash_score_info = cs.info(
      geometry_restraints_manager=grm,
      macro_molecule_selection=macro_mol_sel,
      sites_cart=sites_cart,
      site_labels=site_labels,
      hd_sel=hd_sel)
  return clash_score_info

if __name__ == '__main__':
  run()


