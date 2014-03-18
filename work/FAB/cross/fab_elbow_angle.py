from __future__ import division
from scitbx.linalg import eigensystem
from scitbx.array_family import flex
from libtbx.utils import null_out
from math import acos,pi
from iotbx import pdb
import iotbx.pdb

class fab_elbow_angle(object):

  def __init__(self,
               pdb_hierarchy,
               chain_ID_light='L',
               chain_ID_heavy='H',
               limit_light=107,
               limit_heavy=113):
    '''
    Get elbow angle for Fragment antigen-binding (Fab)

    - Default heavy and light chains IDs are: H : heavy,  L : light
    - Default limit (cutoff) between variable and constant parts
      is residue number 107/113 for light/heavy chains
    - Variable domain is from residue 1 to limit.
      Constant domain form limit+1 to end.
    - Method of calculating angle is based on Stanfield, et al., JMB 2006
    '''
     # create selection strings for the heavy/light var/const part of chains
    self.select_str(
      chain_ID_H=chain_ID_heavy,
      limit_H=limit_heavy,
      chain_ID_L=chain_ID_light,
      limit_L=limit_light)
    # get the hirarchy for and divide using selection strings
    self.pdb_hierarchy = pdb_hierarchy
    self.get_pdb_chains()
    # Get heavy to light reference vector before alignment !!!
    vh_end = self.pdb_var_H.atoms()[-1].xyz
    vl_end = self.pdb_var_L.atoms()[-1].xyz
    mid_H_to_L = self.norm_vec(start=vh_end,end=vl_end)
    # Get transformations objects
    tranformation_const= self.get_transformation(
      fixed_selection=self.pdb_const_H,
      moving_selection=self.pdb_const_L)
    tranformation_var = self.get_transformation(
      fixed_selection=self.pdb_var_H,
      moving_selection=self.pdb_var_L)
    # Get the angle and eigenvalues
    eigen_const = eigensystem.real_symmetric(tranformation_const.r.as_sym_mat3())
    eigen_var = eigensystem.real_symmetric(tranformation_var.r.as_sym_mat3())
    # c : consttant, v : variable
    eigenvectors_c = self.get_eigenvector(eigen_const)
    eigenvectors_v = self.get_eigenvector(eigen_var)
    # test eignevectors pointing in oposite directions
    if eigenvectors_c.dot(eigenvectors_v) > 0:
      eigenvectors_v = - eigenvectors_v
    # Calc Feb elbow angle
    angle = self.get_angle(vec1=eigenvectors_c, vec2=eigenvectors_v)
    # Test if elbow angle larger or smaller than 180
    zaxis = self.cross(eigenvectors_v, eigenvectors_c)
    xaxis = self.cross(eigenvectors_c,zaxis)
    # choose ref axis
    ref_axis = zaxis
    #if abs(mid_H_to_L.dot(xaxis)) > abs(mid_H_to_L.dot(zaxis)):
      #ref_axis = xaxis
    if mid_H_to_L.dot(ref_axis) < 0:
      angle = 360 - angle
    self.fab_elbow_angle = angle

  def norm_vec(self,start,end):
    '''retruns normalized vector that starts at "stat" and ends at "end"'''
    x = flex.double(end) - flex.double(start)
    return x/x.norm()

  def cross(self,a,b):
    '''(array,array) -> array
    returns a normalized cross product vector'''
    a1,a2,a3 = a
    b1,b2,b3 = b
    x = flex.double([a2*b3-a3*b2,a3*b1-a1*b3,a1*b2-a2*b1])
    return x/x.norm()

  def get_angle(self,vec1,vec2,larger=True):
    '''retrun the larger angle between vec1 and vec2'''
    if vec1 and vec1:
      angle_cos = vec1.dot(vec2)
      acos_angle_cos = acos(angle_cos)
      assert acos_angle_cos != 0
      angle = 180/pi*acos_angle_cos
    else:
      angle = 0
    if (angle < 90) and larger: angle = 180 - angle
    if (angle > 90) and not larger: angle = 180 - angle
    return angle

  def get_eigenvector(self,eigen):
    '''
    Get the eigen vector for eigen value 1 and normalize it
    '''
    v = eigen.vectors()
    e = eigen.values()
    indx = None
    # select eigenvector that corespondes to a real egienvalue == 1
    for i,x in enumerate(e):
      if not isinstance(x,complex):
        if abs(1-x)<1e-6:
          indx = i
          break
    # make sure we have egienvalue == 1
    assert not indx
    eigenvector = v[indx:indx+3]
    # normalize
    eigenvector = eigenvector / eigenvector.dot(eigenvector)
    if e.all_eq(flex.double([1,1,1])):
      eigenvector = None
    return eigenvector

  def get_pdb_chains(self):
    '''Create seperate pdb hierarchy for each on the chains we want to align'''
    ph = self.pdb_hierarchy
    # test selection
    test = ph.atom_selection_cache().selection
    #
    self.pdb_var_H = ph.select(test(self.select_var_str_H))
    self.pdb_const_H = ph.select(test(self.select_const_str_H))
    self.pdb_var_L = ph.select(test(self.select_var_str_L))
    self.pdb_const_L = ph.select(test(self.select_const_str_L))

  def get_transformation(self,fixed_selection,moving_selection):
    from phenix.command_line import superpose_pdbs
    '''
    Align the moving pdb hierarchy on to the fixed one.
    Provides an object with rotation and translation info
    '''
    params = superpose_pdbs.master_params.extract()
    x = superpose_pdbs.manager(
      params,
      log=null_out(),
      write_output=False,
      save_lsq_fit_obj=True,
      pdb_hierarchy_fixed=fixed_selection,
      pdb_hierarchy_moving=moving_selection)
    return x.lsq_fit_obj

  def select_str(self,chain_ID_H,limit_H,chain_ID_L,limit_L):
    '''create selection strings for the heavy and light chains
    seperating the vairable and constant parts of the chains'''
    s1 = 'pepnames and (name ca or name n or name c) and altloc " "'
    s2 = 'chain {0} and resseq {1}:{2} and {3}'
    self.select_var_str_H = s2.format(chain_ID_H,1,limit_H,s1)
    self.select_const_str_H = s2.format(chain_ID_H,limit_H+1,'end',s1)
    self.select_var_str_L = s2.format(chain_ID_L,1,limit_L,s1)
    self.select_const_str_L = s2.format(chain_ID_L,limit_L+1,'end',s1)
