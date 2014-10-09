from __future__ import division
import os

from iotbx.pdb.multimer_reconstruction import multimer
import mmtbx.ncs.ncs_utils as nu
import iotbx.ncs
import iotbx.pdb


# ncs_1_copy="""\
# MTRIX1   1  1.000000  0.000000  0.000000        0.00000    1
# MTRIX2   1  0.000000  1.000000  0.000000        0.00000    1
# MTRIX3   1  0.000000  0.000000  1.000000        0.00000    1
# MTRIX1   2  0.496534 -0.643583  0.582456        0.00000
# MTRIX2   2  0.867956  0.376083 -0.324360        0.00000
# MTRIX3   2 -0.010295  0.666605  0.745340        0.00000
# ATOM      1  N   THR A   1      11.797  12.521   4.849  1.00 10.00           N
# ATOM      2  CA  THR A   1      11.607  11.311   4.057  1.00 10.00           C
# ATOM      3  C   THR A   1      11.119  10.155   4.925  1.00 10.00           C
# ATOM      4  O   THR A   1      11.603   9.957   6.039  1.00 10.00           O
# ATOM      5  CB  THR A   1      12.906  10.893   3.344  1.00 10.00           C
# ATOM      6  OG1 THR A   1      13.340  11.949   2.478  1.00 10.00           O
# ATOM      7  CG2 THR A   1      12.683   9.631   2.525  1.00 10.00
# TER
# """

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
ATOM      1  N   THR A   1      11.797  12.521   4.849  1.00 10.00           N
ATOM      2  CA  THR A   1      11.607  11.311   4.057  1.00 10.00           C
ATOM      3  C   THR A   1      11.119  10.155   4.925  1.00 10.00           C
ATOM      4  O   THR A   1      11.603   9.957   6.039  1.00 10.00           O
ATOM      5  CB  THR A   1      12.906  10.893   3.344  1.00 10.00           C
ATOM      6  OG1 THR A   1      13.340  11.949   2.478  1.00 10.00           O
ATOM      7  CG2 THR A   1      12.683   9.631   2.525  1.00 10.00
TER
"""


def run():
  os.chdir(r'C:\Phenix\Dev\Work\work\NCS\junk')

  pdb_inp = iotbx.pdb.input(source_info=None, lines=ncs_1_copy)
  # pdb_obj = iotbx.pdb.hierarchy.input(pdb_string=ncs_1_copy)
  # mtrix_object = pdb_inp.process_mtrix_records()
  ph = pdb_inp.construct_hierarchy()
  xrs = pdb_inp.xray_structure_simple()
  # xrs_one_ncs = xrs.orthorhombic_unit_cell_around_centered_scatterers(
  #   buffer_size=8)
  # ph.adopt_xray_structure(xrs_one_ncs)
  xrs_one_ncs = xrs
  xrs_shaken = xrs_one_ncs.deep_copy_scatterers()
  #
  m = multimer(
    pdb_str=ncs_1_copy,
    round_coordinates=False,
    reconstruction_type='cau',error_handle=True,eps=1e-2)
  xrs_asu = m.assembled_multimer.extract_xray_structure(
    crystal_symmetry = xrs_one_ncs.crystal_symmetry())
  m.write("answer_asu.pdb")
  #
  transforms_obj = iotbx.ncs.input(
    pdb_string=ncs_1_copy)
  x = nu.concatenate_rot_tran(transforms_obj)
  nrg = transforms_obj.get_ncs_restraints_group_list()
  print 'before shaking:'
  print_for_test(x=x, nrg=nrg)
  x = nu.shake_transformations(
    x = x,
    shake_angles_sigma=0.03,
    shake_translation_sigma=0.3)
  transforms_obj = nu.separate_rot_tran(x=x,transforms_obj=transforms_obj)
  nrg = transforms_obj.get_ncs_restraints_group_list()
  print 'after shaking:'
  print_for_test(x=x, nrg=nrg)
  #
  mtrix_object = transforms_obj.build_MTRIX_object()
  ph.adopt_xray_structure(xrs_shaken)
  of = open("master_ncs_shaken.pdb", "w")
  print >> of, mtrix_object.as_pdb_string()
  print >> of, ph.as_pdb_string(crystal_symmetry=xrs.crystal_symmetry())
  of.close()
  #
  m = multimer(
    file_name="master_ncs_shaken.pdb",
    round_coordinates=False,
    reconstruction_type='cau',error_handle=True,eps=1e-2)
  xrs_asu = m.assembled_multimer.extract_xray_structure(
    crystal_symmetry = xrs_one_ncs.crystal_symmetry())
  m.write("shaken_asu.pdb")

  print 'Done'
  print 'Check results in: \nanswer_asu.pdb \nshaken_asu.pdb '

def print_for_test(x,nrg,n=0):
  s = '{:.6f} ' * 6
  xx = nu.concatenate_rot_tran(ncs_restraints_group_list=nrg)
  a1 =  list(x)
  a2 =  list(xx)
  s2 = '{:.6f} ' * 3
  xxx = (x-xx).min_max_mean().as_tuple()
  s1 = '{}: '.format(n)
  s = s1 + s
  print s.format(*a1)
  print s.format(*a2)
  # print s2.format(*xxx)

if __name__ == '__main__':
  run()
