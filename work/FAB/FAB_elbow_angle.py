from __future__ import division
from phenix.command_line import superpose_pdbs
from iotbx import pdb
from iotbx.pdb import fetch
from libtbx import easy_run
from libtbx.utils import Sorry
from libtbx.utils import null_out
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
    '''Get elbow angle for Fragment antigen-binding (Fab)

    If not specified otherwise, the chains IDs are: H : heavy,  L : light

    I assume that the Constant protion of the protein is from start_to_limit
    and that the Vairiable one is from limit_to_end

    !!!!  need to check if the assumption above is true !!!!


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

    '''

    # Test if FAB

    # Devide to variable and constant part, and get the hirarchy for
    # H : heavy,  L : light
    # start_to_limit : Constant
    # limit_to_end : Variable
    pdb_const_H,pdb_var_H = self.get_pdb_protions(
      pdb_file_name=pdb_file_name,
      chain_ID=chain_ID_heavy,
      limit=limit_heavy)
    pdb_const_L,pdb_var_L = self.get_pdb_protions(
      pdb_file_name=pdb_file_name,
      chain_ID=chain_ID_light,
      limit=limit_light)
    # get ratoation and translation information
    tranformation_const = self.get_transformation(
      pdb_hierarchy_fixed=pdb_const_H,
      pdb_hierarchy_moving=pdb_const_L)
    self.rotation_const = tranformation_const.r
    tranformation_var = self.get_transformation(
      pdb_hierarchy_fixed=pdb_var_H,
      pdb_hierarchy_moving=pdb_var_L)
    #
    self.rotation_const = tranformation_const.r
    self.rotation_var = tranformation_var.r
    # Get the angle
    self.FAB_elbow_angle = self.rotation_const.angle(self.rotation_var,deg=True)

  def get_transformation(self,pdb_hierarchy_fixed,pdb_hierarchy_moving):
    '''
    Create a superpose_pdbs manager object, by alinning the fix and moving chains,
    being compared, and calculating transoformation needed to align the moving hierarchy
    with the fixed on.
    It will disregard non-aligned portion of the two, when producing
    rotation and translation transformations.

    Arguments:
    ----------
    pdb_hierarchy_fixed : pdb_hierarchy of a portion of a protein that will stay fix in place
    pdb_hierarchy_moving

    Retrun:
    -------
    lsq_fit_obj : least-squre-fit object that contians the transformation information
    '''
    params = superpose_pdbs.master_params.extract()
    x = superpose_pdbs.manager(
      params,
      log=null_out(),
      #log=None,
      write_output=False,
      save_lsq_fit_obj=True,
      pdb_hierarchy_fixed=pdb_hierarchy_fixed,
      pdb_hierarchy_moving=pdb_hierarchy_moving)
    return x.lsq_fit_obj

  def get_pdb_protions(self,pdb_file_name,chain_ID,limit,write_to_file=False):
    '''(str,str,int,bool) -> pdb_hierarchy,pdb_hierarchy
    Takes a chain, with chain_ID, from a pdb file with pdb_file_name
    and devides it to two, at the limit position.

    Produces two new hierarchies from the two protions of the protein

    Can writes  two pdb parts into new files, in a pdb format
    str1 with the name pdb_file_name + chain_ID + Var .pdb
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
    fn = self.check_pdb_file_name(pdb_file_name)
    pdb_obj = pdb.hierarchy.input(file_name=fn)
    chains = pdb_obj.hierarchy.models()[0].chains()
    chain_names = [x.id for x in chains]
    # check if correct name is given and find the correct chain
    if chain_ID in chain_names:
      indx = chain_names.index(chain_ID)
      chain_to_split = pdb_obj.hierarchy.models()[0].chains()[indx]
      # create copies for the variable and the constant
      pdb_limit_to_end = self.creat_new_pdb_from_chain(
        chain=chain_to_split,
        new_file_name='1bbd_{}_Const.pdb'.format(chain_ID),
        to_end=True,
        write_to_file=write_to_file)
      pdb_start_to_limit = self.creat_new_pdb_from_chain(
        chain=chain_to_split,
        new_file_name='1bbd_{}_Var.pdb'.format(chain_ID),
        to_end=False,
        write_to_file=write_to_file)
    else:
      raise Sorry('The is not chain with {0} name in {1}'.format(chain_ID,fn))
    return pdb_start_to_limit,pdb_limit_to_end

  def creat_new_pdb_from_chain(self,chain,new_file_name,limit=117,to_end=True,write_to_file=False):
    '''(pdb_object,int,bool,bool) -> pdb_hierarchy

    Creates new pdb hierarchy from a portion of chain

    Arguments:
    ----------
    pdb_object : a chain from a pdb file
    limit : the cutoff between the Variable and the Constant FAB domains
    to_end : if False take portion before limit, if true, the protion after
    write_to_file : if True, will create files such as 1bbd_H_Var.pdb, 1bbd_H_Const.pdb

    Returns:
    --------
    pdb_hierarchy1, pdb_hierarchy2
    pdb_hierarchy1: 0 to limit atoms from  pdb_file_name
    pdb_hierarchy2: limit to end atoms from pdb_file_name

    '''
    chain = chain.detached_copy()
    if to_end:
      i_start = 0
      i_end = limit
    else:
      i_start = limit
      i_end = len(chain.residue_groups())
    # Create a new hierarchy, keep chain ID
    new_chain = pdb.hierarchy.chain()
    new_model = pdb.hierarchy.model()
    new_root = pdb.hierarchy.root()
    new_chain.id = chain.id
    for i in range(i_start,i_end):
      residue_group = chain.residue_groups()[i]
      residue_group_copy = residue_group.detached_copy()
      new_chain.append_residue_group(residue_group_copy)
    new_model.append_chain(new_chain)
    new_root.append_model(new_model)
    if write_to_file:
      new_root.write_pdb_file(file_name=new_file_name)
    return new_root

  def check_pdb_file_name(self,file_name):
    tmp = os.path.basename(file_name)
    tmp = tmp.split('.')
    assert len(tmp[0])==4
    if len(tmp)==2 and tmp[1]=='pdb':
      fn = file_name
    else:
      fn = tmp[0] + '.pdb'
    return fn
