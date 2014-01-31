from __future__ import division
from phenix.command_line import superpose_pdbs
from scitbx.linalg import eigensystem
from scitbx.array_family import flex
from libtbx.utils import null_out
from libtbx.utils import Sorry
from iotbx.pdb import fetch
from math import acos,pi
from iotbx import pdb
import shutil
import tempfile
import os,sys


class FAB_elbow_angle(object):
  def __init__(self,
               pdb_file_name,
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
    self.FAB_elbow_angle : the elbow angle calculated as the dot product of
                           the VL-VH pseudodyade axie and the CL-CH pseudodyade axie
                           The angle always computes between 90 and 180

    Test program at:
    cctbx_project\mmtbx\regression\tst_FAB_elbow_angle.py

    Example:
    --------
    >>>fab = FAB_elbow_angle(
         pdb_file_name='1bbd',
         chain_ID_light='L',
         chain_ID_heavy='H',
         limit_light=114,
         limit_heavy=118)
    >>> print fab.FAB_elbow_angle
    133
    >>>fab = FAB_elbow_angle(pdb_file_name='1bbd')
    >>> print fab.FAB_elbow_angle
    126
    (127 in Stanfield, et al., JMB 2006)

    @author Youval Dar (LBL 2014)
    '''
    # abolute path and test that file exist
    pdb_file_name = self.get_pdb_file_name_and_path(pdb_file_name)
    # Devide to variable and constant part, and get the hirarchy for
    # H : heavy,  L : light
    # start_to_limit : Constant
    # limit_to_end : Variable

    select_var_str_H,select_const_str_H  = select_str(chain_ID=chain_ID_heavy,limit=limit_heavy)
    select_var_str_L,select_const_str_L  = select_str(chain_ID=chain_ID_light,limit=limit_light)

    #pdb_var_H,pdb_const_H = self.get_pdb_protions(
      #pdb_file_name=pdb_file_name,
      #chain_ID=chain_ID_heavy,
      #limit=limit_heavy)
    #pdb_var_L,pdb_const_L = self.get_pdb_protions(
      #pdb_file_name=pdb_file_name,
      #chain_ID=chain_ID_light,
      #limit=limit_light)
    # get ratoation and translation information
    tranformation_const = self.get_transformation(
      fixed_selection=select_const_str_H,
      moving_selection=select_const_str_L)
    self.rotation_const = tranformation_const.r
    tranformation_var = self.get_transformation(
      fixed_selection=select_var_str_H,
      moving_selection=select_var_str_L)
    self.rotation_const = tranformation_const.r
    self.rotation_var = tranformation_var.r
    # Get the angle and eigenvalues
    eigen_const = eigensystem.real_symmetric(tranformation_const.r.as_sym_mat3())
    eigen_var = eigensystem.real_symmetric(tranformation_var.r.as_sym_mat3())
    # get 1bbd
    currnet_dir = os.getcwd()
    tempdir = tempfile.mkdtemp('tempdir')
    os.chdir(tempdir)
    fetch.get_pdb ('1bbd',data_type='pdb',mirror='rcsb',log=null_out())
    # align constant heavy of tested protein with a reference protein 1ddb
    test_var_H,test_const_H = self.get_pdb_protions(
      pdb_file_name=pdb_file_name,
      chain_ID=chain_ID_heavy,
      limit=limit_heavy)
    ref_var_H,ref_const_H = self.get_pdb_protions(
      pdb_file_name='1bbd.pdb',
      chain_ID=chain_ID_heavy,
      limit=limit_heavy)
    # clean temp folder
    os.chdir(currnet_dir)
    shutil.rmtree(tempdir)
    # get rotation and translation
    ref_tranformation = self.get_transformation(
      fixed_selection=ref_const_H,
      moving_selection=test_const_H)
    self.ref_rotation = ref_tranformation.r
    self.ref_translation = ref_tranformation.t
    # Apply transformation on the variable portion of the tesed protein
    new_sites = ref_tranformation.r.elems*test_var_H.atoms().extract_xyz() + ref_tranformation.t
    test_var_H.atoms().set_xyz(new_sites)
    # get rotation and translation
    ref_tranformation_var = self.get_transformation(
      fixed_selection=ref_var_H,
      moving_selection=test_var_H)
    self.ref_rotation_var = ref_tranformation_var.r
    # Get the angle and eigenvalues for reference protein
    eigen_ref = eigensystem.real_symmetric(ref_tranformation_var.r.as_sym_mat3())
    eigen_vectors_const = self.get_eigenvector(eigen_const)
    eigen_vectors_var = self.get_eigenvector(eigen_var)
    eigen_vectors_var_ref = self.get_eigenvector(eigen_ref)
    angle = self.get_angle(vec1=eigen_vectors_const, vec2=eigen_vectors_var)
    ref_angle = self.get_angle(vec1=eigen_vectors_var_ref, vec2=eigen_vectors_var,larger=False)
    # Resolve ambiguity with angle
    if ref_angle > 90: ref_angle = 180 - ref_angle
    if angle + ref_angle > 180:
      # Choose angle smaller than 180
      angle = 360 - angle
    self.FAB_elbow_angle = angle

  def get_segment_info(self,pdb_obj,segment_id):
    '''(pdb_hierarchy,str) -> str,int
    Return information on first and last residues in pdb_obj
    as well as the number of atoms in it'''
    segment_start_end = '{0} from {1} to {2}'
    n_atoms_in_segment = len(pdb_obj.atoms())
    s_start = pdb_obj.atoms()[0].pdb_label_columns()
    s_end = pdb_obj.atoms()[-1].pdb_label_columns()
    segment_start_end = segment_start_end.format(segment_id,s_start,s_end)
    return segment_start_end,n_atoms_in_segment

  def get_angle(self,vec1,vec2,larger=True):
    '''retrun the larger angle between vvec1 and vec2'''
    if vec1 and vec1:
      angle_cos = vec1.dot(vec2)
      angle = 180/pi*acos(angle_cos)
    else:
      angle = 0
    if (angle < 90) and larger: angle = 180 - angle
    if (angle > 90) and not larger: angle = 180 - angle
    return angle

  def cross_prod(self,a,b):
    '''(array,array) -> array'''
    a1,a2,a3 = a
    b1,b2,b3 = b
    return flex.double([a2*b3-a3*b2,a3*b1-a1*b3,a1*b2-a2*b1])

  def get_eigenvector(self,eigen):
    '''
    Get the eigen vector for eigen value 1
    and normalize it
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

  def get_transformation(self,fixed_selection,moving_selection):
    '''
    Create a superpose_pdbs manager object, by alinning the fix and moving chains,
    being compared, and calculating transoformation needed to align the moving hierarchy
    with the fixed on.
    It will disregard non-aligned portion of the two, when producing
    rotation and translation transformations.

    Arguments:
    ----------
    fixed_selection : pdb_hierarchy of a portion of a protein that will stay fix in place
    moving_selection

    Retrun:
    -------
    lsq_fit_obj : least-squre-fit object that contians the transformation information
    '''
    params = superpose_pdbs.master_params.extract()
    params.selection_default_fixed = pdb_hierarchy_fixed
    params.selection_default_moving = moving_selection
    x = superpose_pdbs.manager(
      params,
      log=null_out(),
      write_output=False,
      save_lsq_fit_obj=True)
    return x.lsq_fit_obj

  def get_pdb_protions(self,pdb_file_name,chain_ID,limit,write_to_file=False):
    '''(str,str,int,bool) -> pdb_hierarchy,pdb_hierarchy
    Takes a chain, with chain_ID, from a pdb file with pdb_file_name
    and devides it to two, at the limit position.

    Produces two new hierarchies from the two protions of the protein

    Can writes  two pdb parts into new files, in a pdb format
    (Start to Limit -> Variable)
    str1 with the name pdb_file_name + chain_ID + Var .pdb
    (Limit to End -> Constant)
    str2 with the name pdb_file_name + chain_ID + Const .pdb

    Example:
    >>>pdb_start_to_limit, pdb_limit_to_end = get_pdb_protions(pdb_file_name='1bbd',chain_ID='H',limit=113)


    Arguments:
    ----------
    pdb_file_name : pdb file name
    chain_ID : PDB chain ID
    limit : the cutoff between the Variable and the Constant FAB domains
    write_to_file : if True, will create files such as 1bbd_H_Var.pdb, 1bbd_H_Const.pdb

    Returns:
    --------
    pdb_hierarchy1, pdb_hierarchy2
    pdb_hierarchy1: 0 to limit atoms from  pdb_file_name
    pdb_hierarchy2: limit to end atoms from pdb_file_name
    '''
    pdb_hierarchy = pdb.hierarchy.input(file_name=pdb_file_name).hierarchy

    ss = 'pepnames and (name ca or name n or name c) and altloc " "'

    select_var_str = 'chain {0} and resseq {1}:{2} and {3}'.format(chain_ID,1,limit,ss)
    select_const_str = 'chain {0} and resseq {1}:end and {2}'.format(chain_ID,limit+1,ss)
    selection_var = pdb_hierarchy.atom_selection_cache().selection(select_var_str)
    selection_const = pdb_hierarchy.atom_selection_cache().selection(select_const_str)

    pdb_var = pdb_hierarchy.select(selection_var)
    pdb_const = pdb_hierarchy.select(selection_const)

    if write_to_file:
      fn = os.path.splitext(os.path.basename(pdb_file_name))[0]
      var_file_name='{0}_{1}_Var.pdb'.format(fn,chain_ID)
      pdb_const='{0}_{1}_Const.pdb'.format(fn,chain_ID)
      pdb_var.write_pdb_file(file_name=var_file_name)
      pdb_const.write_pdb_file(file_name=pdb_const)

    return pdb_var,pdb_const

  def select_str(chain_ID,limit):
    ss = 'pepnames and (name ca or name n or name c) and altloc " "'
    select_var_str = 'chain {0} and resseq {1}:{2} and {3}'.format(chain_ID,1,limit,ss)
    select_const_str = 'chain {0} and resseq {1}:end and {2}'.format(chain_ID,limit+1,ss)
    return select_var_str,select_const_str


  def get_pdb_file_name_and_path(self,file_name):
    tmp = os.path.basename(file_name)
    tmp = tmp.split('.')
    assert len(tmp[0])==4
    if len(tmp)==2 and tmp[1]=='pdb':
      fn = file_name
    else:
      fn = tmp[0] + '.pdb'
    assert os.path.isfile(fn)
    return os.path.abspath(fn)
