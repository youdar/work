from __future__ import division
from scitbx.linalg import eigensystem
from scitbx.array_family import flex
from libtbx.utils import null_out
from math import acos,pi
from scitbx import matrix
from iotbx import pdb
import os


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
    - Variable domain si from residue 1 to limit.
      Constant domain form limit+1 to end.
    - Method of calculating angle is based on Stanfield, et al., JMB 2006

    Argument:
    ---------
    pdb_file_name : 4 characters string, a PDB name
    chain_ID_heavy : The heavy protion of the protein, chain ID
    chain_ID_light : The light protion of the protein, chain ID
    limit_heavy : the number of the cutoff residue, between
                  the variable and constant portions in the heavy chian
    limit_light : the number of the cutoff residue, between
                  the variable and constant portions in the light chian

    Main attributes:
    ----------------
    self.fab_elbow_angle : the elbow angle calculated as the dot product of
                           the VL-VH pseudodyade axie and the CL-CH
                           pseudodyade axie

    Test program at:
    cctbx_project\mmtbx\regression\tst_fab_elbow_angle.py

    Example:
    --------
    >>>fab = fab_elbow_angle(
         pdb_file_name='1bbd',
         chain_ID_light='L',
         chain_ID_heavy='H',
         limit_light=114,
         limit_heavy=118)
    >>> print fab.fab_elbow_angle
    133
    >>>fab = fab_elbow_angle(pdb_file_name='1bbd')
    >>> print fab.fab_elbow_angle
    126 (127 in Stanfield, et al., JMB 2006)

    @author Youval Dar (LBL 2014)
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
    #mid_H_to_L = self.H_to_L_vec()

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

    # c: res 171    v: res 44
    c = 20*eigenvectors_c + flex.double(self.pdb_const_H.atoms()[135].xyz)
    v = 20*eigenvectors_v + flex.double(self.pdb_var_H.atoms()[130].xyz)
    r = 20*mid_H_to_L + flex.double(self.pdb_var_H.atoms()[-1].xyz)
    rs = flex.double(self.pdb_var_H.atoms()[-1].xyz)
    re = flex.double(self.pdb_var_L.atoms()[-1].xyz)
    print
    print('c')
    print list(flex.double(self.pdb_const_H.atoms()[135].xyz))
    print list(c)
    print 'v'
    print list(flex.double(self.pdb_var_H.atoms()[130].xyz))
    print list(v)
    print 'r'
    print list(flex.double(self.pdb_var_H.atoms()[-1].xyz))
    print list(r)
    print 'f'
    print self.pdb_var_H.atoms()[-1].id_str()
    print list(rs)
    print self.pdb_var_L.atoms()[-1].id_str()
    print list(re)



    #


    # test eignevectors pointing in oposite directions
    if eigenvectors_c.dot(eigenvectors_v) > 0:
      print 'reversing direction of variable rotation eigenvector!!!!'
      eigenvectors_v = - eigenvectors_v
    # Calc Feb elbow angle
    angle = self.get_angle(vec1=eigenvectors_c, vec2=eigenvectors_v)
    # Test if elbow angle larger or smaller than 180
    zaxis = self.cross(eigenvectors_v, eigenvectors_c)
    xaxis = self.cross(eigenvectors_c,zaxis)

    x = 20*xaxis + flex.double(self.pdb_const_H.atoms()[135].xyz)
    print 'x'
    print list(flex.double(self.pdb_const_H.atoms()[135].xyz))
    print list(x)

    print mid_H_to_L.dot(xaxis)
    print mid_H_to_L.dot(zaxis)

    #m = matrix.sqr(list(xaxis)+list(eigenvectors_c)+list(zaxis))
    #A = m.transpose()
    #Ainv = A.inverse(0
    # choose ref axis
    ref_axis = zaxis
    #if abs(mid_H_to_L.dot(xaxis)) > abs(mid_H_to_L.dot(zaxis)):
      #ref_axis = xaxis
    if mid_H_to_L.dot(ref_axis) < 0:
      angle = 360 - angle
    self.fab_elbow_angle = angle

  def H_to_L_vec(self):
    ''' get the vector from the center of coordinates of the heavy chain
    to the center of coordinates of the light chain'''
    H = flex.double([0,0,0])
    L = flex.double([0,0,0])
    for x in self.pdb_const_H.atoms(): H += flex.double(x.xyz)
    for x in self.pdb_var_H.atoms(): H += flex.double(x.xyz)
    for x in self.pdb_const_L.atoms(): L += flex.double(x.xyz)
    for x in self.pdb_var_L.atoms(): L += flex.double(x.xyz)
    H = H/(len(self.pdb_const_H.atoms())+len(self.pdb_var_H.atoms()))
    L = L/(len(self.pdb_const_L.atoms())+len(self.pdb_var_L.atoms()))
    return self.norm_vec(start=H,end=L)

  def norm_vec(self,start,end):
    ''' retruns normalized vector that starts at "stat" and ends at "end"'''
    x = flex.double(end) - flex.double(start)
    return x/x.norm()


  def get_angle(self,vec1,vec2,larger=True):
    '''retrun the larger angle between vec1 and vec2'''
    if vec1 and vec1:
      angle_cos = vec1.dot(vec2)
      angle = 180/pi*acos(angle_cos)
    else:
      angle = 0
    if (angle < 90) and larger: angle = 180 - angle
    if (angle > 90) and not larger: angle = 180 - angle
    return angle

  def cross(self,a,b):
    '''(array,array) -> array
    returns a normalized cross product vector'''
    a1,a2,a3 = a
    b1,b2,b3 = b
    x = flex.double([a2*b3-a3*b2,a3*b1-a1*b3,a1*b2-a2*b1])
    return x/x.norm()

  def get_eigenvector(self,eigen):
    '''
    Get the eigen vector for eigen value 1
    and normalize it
    '''
    v = eigen.vectors()
    e = eigen.values()
    indx = None
    l = 0
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

    Arguments:
    ----------
    fixed_selection, moving_selection : pdb_hierarchy

    Retrun:
    -------
    lsq_fit_obj : least-squre-fit object that contians the
                  transformation information
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
    sel_str = lambda ID,i_s,i_e,s: s2.format(ID,i_s,i_e,s)
    self.select_var_str_H = sel_str(chain_ID_H,1,limit_H,s1)
    self.select_const_str_H = sel_str(chain_ID_H,limit_H+1,'end',s1)
    self.select_var_str_L = sel_str(chain_ID_L,1,limit_L,s1)
    self.select_const_str_L = sel_str(chain_ID_L,limit_L+1,'end',s1)

