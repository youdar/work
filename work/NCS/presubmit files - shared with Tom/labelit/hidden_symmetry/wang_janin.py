import math

from iotbx import pdb
from libtbx.test_utils import approx_equal
import scitbx.math
from scitbx import matrix
from phenix.command_line.simple_ncs_from_pdb import simple_ncs_from_pdb

def compare_axes(ncs_pdb_file=None,crystal_symmetry=None,cb_op=None,eps_m=0.05):
  """
  =============================================================================
  Function calculates the angular difference between NCS and crystallographic
  symmetry axes as described in,

    Wang & Janin, Acta Cryst. D49, 505-512.

  Arguments:
    ncs_pdb_file - the PDB file name for input into phenix.simple_ncs_from_pdb
    crystal_symmetry - the crystal symmetry for comparison (crystal.symmetry)
    cb_op - change of basis that transforms the setting of crystal_symmetry
            to the setting in ncs_pdb_file
    eps_m - tolerance for accepting a rotation matrix as the identity matrix

  Returns:
    A dictionary listing the axis and magnitude for each NCS operator and its
    angular difference (degrees) from a crystallographic symmetry axis.

  Notes:
    The setting for ncs_pdb_file and crystal_symmetry should be the same (by
    using cb_op).
    The key for each entry is also the angular difference for sorting purposes.
  -----------------------------------------------------------------------------
  """

  # find ncs in file
  ncs = simple_ncs_from_pdb(args=[ncs_pdb_file],quiet=True,suppress_print=True)
  if (len(ncs.ncs_object.ncs_groups()) == 0):
    return None
  matrices = []
  I = (1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0)
  for group in ncs.ncs_object.ncs_groups():
    for m in group.rota_matrices():
      if (approx_equal(m,I,eps=eps_m,out=None) is not True):
        matrices.append(m)

  # compare symmetry matrices with ncs rotation matrices
  uc = crystal_symmetry.unit_cell().change_basis(cb_op)
  om = matrix.sqr(uc.orthogonalization_matrix())
  fm = matrix.sqr(uc.fractionalization_matrix())
  answer_dict = {}
  used_matrix = []

  for ncs_matrix in matrices:
    for sym_op in crystal_symmetry.space_group():
      # convert rotation matrix from fractional space to Cartesian space
      sym_op_matrix = (om*matrix.sqr(cb_op.apply(sym_op).r().as_double())*\
                       fm).as_float()

      # compare elements in NCS and sym_op matrices
      if (approx_equal(ncs_matrix,sym_op_matrix,eps=eps_m,out=None)):

        n = scitbx.math.r3_rotation_axis_and_angle_from_matrix(ncs_matrix)
        ncs_axis = matrix.row(n.axis).normalize()
        ncs_angle = n.angle(deg=True)

        s = scitbx.math.r3_rotation_axis_and_angle_from_matrix(sym_op_matrix)
        sym_op_axis = matrix.row(s.axis).normalize()
        sym_op_angle = s.angle(deg=True)

        dp = math.fabs(sym_op_axis.dot(ncs_axis))
        axis_separation = math.degrees(math.acos(dp))
        answer_dict[axis_separation]= [
          '(%6.4f %6.4f %6.4f)'%(ncs_axis[0],ncs_axis[1],ncs_axis[2]),
          ncs_angle,axis_separation]

  return answer_dict

def find_pseudorotations(ncs_pdb_file=None,crystal_symmetry=None,cb_op=None,
                         eps=0.05):
  """
  =============================================================================
  Function finds NCS rotation operators that are pseudosymmetric

  Arguments:
    ncs_pdb_file - the PDB file name for input into phenix.simple_ncs_from_pdb
    crystal_symmetry - the crystal symmetry for comparison (crystal.symmetry)
    cb_op - change of basis that transforms the setting of crystal_symmetry
            to the setting in ncs_pdb_file
    eps - tolerance for a matrix being the identity (per element tolerance)

  Returns:
    A list containing the pseudosymmetric NCS rotation matrices and its
    idealized crystallographic match

  Notes:
    The setting for ncs_pdb_file and crystal_symmetry should be the same (by
    using cb_op).
    The symmetry in ncs_pdb_file should be a subgroup of crystal_symmetry
    The crystal_symmetry symmetry operators missing from ncs_pdb_file are used
    for determining pseudosymmetry (the NCS operators with the smallest
    difference from the missing symmetry operators are returned)
  -----------------------------------------------------------------------------
  """

  # find ncs in file
  ncs = simple_ncs_from_pdb(args=[ncs_pdb_file],
                            quiet=True,suppress_print=True)
  if (len(ncs.ncs_object.ncs_groups()) == 0):
    return None
  matrices = []
  I = (1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0)
  for group in ncs.ncs_object.ncs_groups():
    for m in group.rota_matrices():
      if (approx_equal(m,I,eps=eps,out=None) is not True):
        matrices.append(m)

  # find missing symmetry operators in ncs_pdb_file
  missing_symops = []
  original_symmetry = pdb.input(file_name=ncs_pdb_file).crystal_symmetry()
  for symop0 in crystal_symmetry.space_group():
    new_symop = cb_op.apply(symop0)
    count = 0;
    for symop1 in original_symmetry.space_group():
      if (new_symop.r() == symop1.r()):
        count = count + 1
    if (count < 1):
      missing_symops.append(new_symop)

  # compare symmetry matrices with ncs rotation matrices
  uc = crystal_symmetry.unit_cell().change_basis(cb_op)
  om = matrix.sqr(uc.orthogonalization_matrix())
  fm = matrix.sqr(uc.fractionalization_matrix())

  result = []
  for symop in missing_symops:
    # convert rotation matrix from fractional space to Cartesian space
    symop_matrix = (om*matrix.sqr(symop.r().as_double())*fm).as_float()

    # find ncs operator closest to symmetry operator (rotation only)
    ncs_rmsd = {}
    for ncs_matrix in matrices:
      rmsd = matrix_rmsd(m1=ncs_matrix,m2=symop_matrix)
      ncs_rmsd[rmsd] = {'ncs':ncs_matrix,'cs':symop_matrix,'rmsd':rmsd}

    # pick ncs operator with lowest rmsd
    sorted_keys = sorted(ncs_rmsd.keys())
    result.append(ncs_rmsd[sorted_keys[0]])

  return result

def matrix_rmsd(m1=None,m2=None):
  """
  =============================================================================
  Function returns the rmsd between two matrices by calculating the difference
  for each element in the matrix,

    rmsd = sqrt((1/n)*sum((e1 - e2)^2))

  where n is the number of elements, the sum is over all the elements in a
  matrix, e1 is an element from m1, and e2 is an element from m2.

  Arguments:
    m1 - first matrix (tuple)
    m2 - second matrix (tuple)

  Returns:
    The rmsd between two matrices (float)

  Notes:
    m1 and m2 can be any 1d iterable type that returns its length via len()
    command
  -----------------------------------------------------------------------------
  """
  assert (len(m1) == len(m2))
  n = len(m1)
  rmsd = 0.0
  for i in xrange(n):
    d = m1[i] - m2[i]
    rmsd = rmsd + d*d
  rmsd = math.sqrt((1.0/float(n))*rmsd)
  return rmsd
