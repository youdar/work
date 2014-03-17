from __future__ import division
from iotbx.pdb.multimer_reconstruction import multimer
import iotbx.pdb
from scitbx.array_family import flex

ncs_1_copy="""\
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

def run():
  # 1 NCS copy: starting template to generate whole asu; place into P1 box
  pdb_inp = iotbx.pdb.input(source_info=None, lines=ncs_1_copy)
  mtrix_object = pdb_inp.process_mtrix_records()
  ph = pdb_inp.construct_hierarchy()
  xrs = pdb_inp.xray_structure_simple()
  f = open("one_ncs_in_asu.pdb", "w")
  print >> f, mtrix_object.format_MTRIX_pdb_string()
  print >> f, ph.as_pdb_string(crystal_symmetry=xrs.crystal_symmetry())
  f.close()
  # get the coordinates of the single NCS
  ncs_xyz = pdb_inp.atoms().extract_xyz().as_double()
  # 1 NCS copy -> full asu (expand NCS). This is the answer-strucure
  m = multimer("one_ncs_in_asu.pdb",'cau',error_handle=True,eps=1e-2)
  # m.assembled_multimer is a iotbx_pdb_hierarchy_ext.root
  m_ph = m.assembled_multimer
  xrs_asu = m_ph.extract_xray_structure()

  xyz_m = m.sites_cart()
  xyz_xrs = xrs_asu.sites_cart()

  frac = xrs_asu.unit_cell().fractionalization_matrix()
  orth = xrs_asu.unit_cell().orthogonalization_matrix()

  rotations = m.rotation_matrices
  translations =  m.translation_vectors
  new_x = list(ncs_xyz)
  x = flex.vec3_double(ncs_xyz)
  for r,t in zip(rotations,translations):
    tmp_x = r.elems*x + t
    new_x += list(tmp_x.as_double())
  xyz_direct = flex.vec3_double(flex.double(new_x))

  xrs_asu.set_sites_cart(xyz_m)

  d1 = flex.sqrt((xyz_m - xyz_xrs).dot()).as_double().min_max_mean().as_tuple()
  print 'xyz_m - xyz_xrs'
  print d1
  print '-'*50

  d2 = flex.sqrt((xyz_m - xyz_direct).dot()).as_double().min_max_mean().as_tuple()
  print 'xyz_m - xyz_direct'
  print d2
  print '-'*50

  d3 = flex.sqrt((xyz_xrs - xyz_direct).dot()).as_double().min_max_mean().as_tuple()
  print 'xyz_xrs - xyz_direct'
  print d3
  print '-'*50



if __name__=='__main__':
  run()