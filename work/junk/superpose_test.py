from __future__ import division
import mmtbx.utils.ncs_utils as nu
from scitbx import matrix
from scitbx.math import superpose
from scitbx.array_family import flex
from iotbx import pdb




def test_superpos_pdb():
  """
  verify creation of transformations using superpose_pdb
  """
  # read file and create pdb object
  pdb_obj = pdb.hierarchy.input(pdb_string=pdb_test_data1)
  temp = pdb_obj.hierarchy.atom_selection_cache()

  # set selection strings
  s1 = 'chain A'
  s2 = 'chain C'

  pdb_master_ncs = pdb_obj.hierarchy.select(temp.selection(s1))
  pdb_ncs_copy = pdb_obj.hierarchy.select(temp.selection(s2))

  r = matrix.sqr([0.309017, -0.809017, 0.500000,
                  0.809017,  0.500000, 0.309017,
                  -0.500000, 0.309017, 0.809017])
  t = matrix.col([0,0,0])

  lsq_fit_obj = superpose.least_squares_fit(
    reference_sites = pdb_ncs_copy.atoms().extract_xyz(),
    other_sites     = pdb_master_ncs.atoms().extract_xyz())

  rotation_matrices = lsq_fit_obj.r
  translation_vectors = lsq_fit_obj.t

  x = flex.vec3_double([(28.392,67.262,97.682)])
  y1 = r.elems*x + t.elems
  y2 = rotation_matrices.elems*x + translation_vectors.elems

  print '---------- rotation matrix elems ----------'
  outstr = '%.6f   ' * 9
  print outstr % (r.elems)
  print outstr % (rotation_matrices.elems)
  print
  print '---------- rotation angles -------------'
  print nu.rotation_to_angles(r).as_numpy_array()
  print nu.rotation_to_angles(rotation_matrices).as_numpy_array()
  print
  print '----------  translation ---------- '
  outstr = ' %.6f   ' * 3
  print outstr % t.elems
  print outstr % translation_vectors.elems
  print
  print '---------- coordinates after transform ------'
  print y1.as_double().as_numpy_array()
  print y2.as_double().as_numpy_array()


pdb_test_data1="""\
ATOM    749  O   UNK A  90      28.392  67.262  97.682  1.00  0.00           O
ATOM    750  N   UNK A  91      30.420  66.924  98.358  1.00  0.00           N
TER
ATOM   1495  N   UNK B  67      33.124   2.704 114.920  1.00  0.00           N
TER
ATOM    749  O   UNK C  90       3.199  86.786  85.616  1.00  0.00           O
ATOM    750  N   UNK C  91       4.437  88.467  85.044  1.00  0.00           N
TER
ATOM   1495  N   UNK D  67      65.508  63.662  77.246  1.00  0.00           N
"""

if __name__ == "__main__":
  test_superpos_pdb()







