from __future__ import division
from iotbx.pdb.multimer_reconstruction import multimer
from iotbx import pdb
import iotbx.pdb

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
ATOM      1  N   THR A   1       9.670  10.289  11.135  1.00 20.00           N
ATOM      2  CA  THR A   1       9.559   8.931  10.615  1.00 20.00           C
ATOM      3  C   THR A   1       9.634   7.903  11.739  1.00 20.00           C
ATOM      4  O   THR A   1      10.449   8.027  12.653  1.00 20.00           O
ATOM      5  CB  THR A   1      10.660   8.630   9.582  1.00 20.00           C
ATOM      6  OG1 THR A   1      10.560   9.552   8.490  1.00 20.00           O
ATOM      7  CG2 THR A   1      10.523   7.209   9.055  1.00 20.00           C
TER
"""

def run():
  # 1 NCS copy: starting template to generate whole asu; place into P1 box
  pdb_inp = iotbx.pdb.input(source_info=None, lines=ncs_1_copy)
  mtrix_object = pdb_inp.process_mtrix_records()
  ph = pdb_inp.construct_hierarchy()
  xrs = pdb_inp.xray_structure_simple()
  xrs_one_ncs = xrs.orthorhombic_unit_cell_around_centered_scatterers(
    buffer_size=8)
  ph.adopt_xray_structure(xrs)
  of = open("one_ncs_in_asu.pdb", "w")
  print >> of, mtrix_object.format_MTRIX_pdb_string()
  print >> of, ph.as_pdb_string(crystal_symmetry=xrs.crystal_symmetry())
  of.close()
  # 1 NCS copy -> full asu (expand NCS). This is the answer-strucure
  m = multimer("one_ncs_in_asu.pdb",'cau',error_handle=True,eps=1e-2)
  assert m.number_of_transforms == 3
  xrs_asu = m.assembled_multimer.extract_xray_structure(
    crystal_symmetry = xrs.crystal_symmetry())
  m.write("full_asu.pdb")
  # Generate Fobs from answer structure
  f_obs = abs(xrs_asu.structure_factors(d_min=2, algorithm="direct").f_calc())
  r_free_flags = f_obs.generate_r_free_flags()
  mtz_dataset = f_obs.as_mtz_dataset(column_root_label="F-obs")
  mtz_dataset.add_miller_array(
    miller_array=r_free_flags,
    column_root_label="R-free-flags")
  mtz_object = mtz_dataset.mtz_object()
  mtz_object.write(file_name = "data.mtz")
  # Shake structure - subject to refinement input
  xrs_shaken = xrs_one_ncs.deep_copy_scatterers()
  xrs_shaken.shake_sites_in_place(mean_distance=0.3)
  ph.adopt_xray_structure(xrs_shaken)
  of = open("one_ncs_in_asu_shaken.pdb", "w")
  print >> of, mtrix_object.format_MTRIX_pdb_string()
  print >> of, ph.as_pdb_string(crystal_symmetry=xrs.crystal_symmetry())
  # Refinement

if __name__ == "__main__":
  run()
