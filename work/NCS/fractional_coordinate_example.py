from __future__ import division
import iotbx.pdb
from scitbx.array_family import flex

pdb_str="""\
CRYST1   15.000   15.000   15.000  90.00 90.00  90.00 P 1           1
ATOM      1  CB  PHE H   1       7.767   5.853   7.671  1.00 15.00           C
ATOM      2  CG  PHE H   1       6.935   5.032   8.622  1.00 15.00           C
ATOM      3  CD1 PHE H   1       5.918   4.176   8.140  1.00 15.00           C
ATOM      4  CD2 PHE H   1       7.161   5.107  10.012  1.00 15.00           C
ATOM      5  CE1 PHE H   1       5.126   3.395   9.038  1.00 15.00           C
ATOM      6  CE2 PHE H   1       6.382   4.336  10.930  1.00 15.00           C
ATOM      7  CZ  PHE H   1       5.360   3.476  10.439  1.00 15.00           C
ATOM      8  C   PHE H   1       7.956   7.811   6.133  1.00 15.00           C
ATOM      9  O   PHE H   1       8.506   7.237   5.169  1.00 15.00           O
ATOM     10  OXT PHE H   1       8.143   9.010   6.428  1.00 15.00           O
ATOM     11  H1  PHE H   1       6.253   5.840   5.439  1.00 15.00           H
ATOM     12  H2  PHE H   1       5.364   7.253   5.745  1.00 15.00           H
ATOM     13  N   PHE H   1       5.875   6.461   6.183  1.00 15.00           N
ATOM     14  H3  PHE H   1       5.216   5.927   6.782  1.00 15.00           H
ATOM     15  CA  PHE H   1       7.000   7.000   7.000  1.00 15.00           C
"""

def run():
  pdb_inp = iotbx.pdb.input(source_info=None, lines=pdb_str)
  ph = pdb_inp.construct_hierarchy()
  xray_structure = pdb_inp.xray_structure_simple()
  sites_cart = xray_structure.sites_cart()
  sites_cart3 = ph.atoms().extract_xyz()
  sites_frac = xray_structure.sites_frac()
  #
  sites_cart2 = xray_structure.unit_cell().orthogonalization_matrix()*sites_frac
  print flex.sqrt((sites_cart-sites_cart2).dot()).min_max_mean().as_tuple()
  #
  sites_frac2 = xray_structure.unit_cell().fractionalization_matrix()*sites_cart
  print flex.sqrt((sites_frac-sites_frac2).dot()).min_max_mean().as_tuple()
  #
  print flex.sqrt((sites_cart3-sites_cart).dot()).min_max_mean().as_tuple()


if __name__=='__main__':
  run()