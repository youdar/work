from __future__ import division
from iotbx.pdb.multimer_reconstruction import multimer
import iotbx.pdb
from scitbx.array_family import flex

ncs="""\
MTRIX1   1  1.000000  0.000000  0.000000        0.00000    1
MTRIX2   1  0.000000  1.000000  0.000000        0.00000    1
MTRIX3   1  0.000000  0.000000  1.000000        0.00000    1
MTRIX1   2  0.496590 -0.643597  0.582393        0.00000
MTRIX2   2  0.867925  0.376088 -0.324443        0.00000
MTRIX3   2 -0.010221  0.666588  0.745356        0.00000
MTRIX1   3 -0.317946 -0.173437  0.932111        0.00000
MTRIX2   3  0.760735 -0.633422  0.141629        0.00000
MTRIX3   3  0.565855  0.754120  0.333333        0.00000
ATOM      1  N   THR A   1       8.111  11.079  10.645  1.00 20.00           N
ATOM      2  CA  THR A   1       8.000   9.721  10.125  1.00 20.00           C
ATOM      3  C   THR A   1       8.075   8.693  11.249  1.00 20.00           C
ATOM      4  O   THR A   1       8.890   8.817  12.163  1.00 20.00           O
ATOM      5  CB  THR A   1       9.101   9.420   9.092  1.00 20.00           C
ATOM      6  OG1 THR A   1       9.001  10.342   8.000  1.00 20.00           O
ATOM      7  CG2 THR A   1       8.964   8.000   8.565  1.00 20.00           C
TER
"""

asu = """\
CRYST1   17.101   19.080   20.163  90.00  90.00  90.00 P 1
SCALE1      0.058476  0.000000  0.000000        0.00000
SCALE2      0.000000  0.052411  0.000000        0.00000
SCALE3      0.000000  0.000000  0.049596        0.00000
ATOM      1  N   THR A   1       8.111  11.080  10.645  1.00 20.00           N
ATOM      2  CA  THR A   1       8.000   9.722  10.125  1.00 20.00           C
ATOM      3  C   THR A   1       8.075   8.694  11.249  1.00 20.00           C
ATOM      4  O   THR A   1       8.890   8.818  12.163  1.00 20.00           O
ATOM      5  CB  THR A   1       9.101   9.421   9.092  1.00 20.00           C
ATOM      6  OG1 THR A   1       9.001  10.343   8.000  1.00 20.00           O
ATOM      7  CG2 THR A   1       8.964   8.000   8.565  1.00 20.00           C
TER
ATOM      1  N   THRaa   1       3.096   7.753  15.237  1.00 20.00           N
ATOM      2  CA  THRaa   1       3.612   7.315  13.946  1.00 20.00           C
ATOM      3  C   THRaa   1       4.966   6.629  14.097  1.00 20.00           C
ATOM      4  O   THRaa   1       5.823   7.086  14.853  1.00 20.00           O
ATOM      5  CB  THRaa   1       3.751   8.492  12.964  1.00 20.00           C
ATOM      6  OG1 THRaa   1       2.472   9.107  12.765  1.00 20.00           O
ATOM      7  CG2 THRaa   1       4.291   8.010  11.625  1.00 20.00           C
TER
ATOM      1  N   THRab   1       5.422   0.660  16.494  1.00 20.00           N
ATOM      2  CA  THRab   1       5.208   1.362  15.233  1.00 20.00           C
ATOM      3  C   THRab   1       6.410   2.229  14.875  1.00 20.00           C
ATOM      4  O   THRab   1       6.981   2.900  15.735  1.00 20.00           O
ATOM      5  CB  THRab   1       3.947   2.244  15.285  1.00 20.00           C
ATOM      6  OG1 THRab   1       2.801   1.429  15.560  1.00 20.00           O
ATOM      7  CG2 THRab   1       3.746   2.965  13.960  1.00 20.00           C
TER
"""

def run():
  pdb_inp1 = iotbx.pdb.input(source_info=None, lines=asu)
  ph = pdb_inp1.construct_hierarchy()
  xrs1 = pdb_inp1.xray_structure_simple()


  sites_cart1 = ph.atoms().extract_xyz()
  sites_cart2 = xrs1.sites_cart()
  d = flex.sqrt((sites_cart1 - sites_cart2).dot()).as_double().min_max_mean().as_tuple()
  print 'sites cart pdb - sites cart xrs'
  print d
  print '-'*50


  # 1 NCS copy: starting template to generate whole asu; place into P1 box
  pdb_inp2 = iotbx.pdb.input(source_info=None, lines=ncs)
  mtrix_object = pdb_inp2.process_mtrix_records()
  ph2 = pdb_inp2.construct_hierarchy()
  xrs2 = pdb_inp2.xray_structure_simple()
  f = open("one_ncs_in_asu.pdb", "w")
  print >> f, mtrix_object.format_MTRIX_pdb_string()
  print >> f, ph.as_pdb_string(crystal_symmetry=xrs2.crystal_symmetry())
  f.close()
  # 1 NCS copy -> full asu (expand NCS). This is the answer-strucure
  m = multimer("one_ncs_in_asu.pdb",'cau',error_handle=True,eps=1e-2)
  m.write('test_asu.pdb')
  # m.assembled_multimer is a iotbx_pdb_hierarchy_ext.root
  m_ph = m.assembled_multimer
  xrs_asu = m_ph.extract_xray_structure()

  frac = xrs_asu.unit_cell().fractionalization_matrix()
  orth = xrs_asu.unit_cell().orthogonalization_matrix()

  sites_cart_m = m.sites_cart()
  sites_cart = xrs_asu.sites_cart()
  frac_xrs = xrs_asu.sites_frac()
  frac_m = frac*sites_cart_m
  sites_cart_xrs = orth*frac_xrs

  pdb_inp = iotbx.pdb.input(file_name='test_asu.pdb')
  xray_structure = pdb_inp.xray_structure_simple()
  sites_cart2 = xray_structure.sites_cart()

  d1 = flex.sqrt((frac_m - frac_xrs).dot()).as_double().min_max_mean().as_tuple()
  print 'frac_m - frac_xrs'
  print d1
  print '-'*50

  d2 = flex.sqrt((sites_cart_m - sites_cart_xrs).dot()).as_double().min_max_mean().as_tuple()
  print 'sites_cart_m - sites_cart_xrs'
  print d2
  print '-'*50

  d3 = flex.sqrt((sites_cart_m - sites_cart).dot()).as_double().min_max_mean().as_tuple()
  print 'sites_cart_m - sites_cart'
  print d3
  print '-'*50

  d4 = flex.sqrt((sites_cart_m - sites_cart2).dot()).as_double().min_max_mean().as_tuple()
  print 'sites_cart_m - sites_cart2'
  print d4
  print '-'*50

  d5 = flex.sqrt((sites_cart - sites_cart2).dot()).as_double().min_max_mean().as_tuple()
  print 'sites_cart - sites_cart2'
  print d5
  print '-'*50



if __name__=='__main__':
  run()