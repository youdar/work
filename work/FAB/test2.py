from __future__ import division
from phenix.command_line import superpose_pdbs
from iotbx import pdb
from iotbx.pdb import fetch
from libtbx import easy_run
from libtbx.utils import Sorry
import shutil
import tempfile
import os,sys

def run(fn):
  ''' '''
  # Test if FAB

  # Devide to
  pdb_hierarchy_var_H,pdb_hierarchy_const_H = get_pdb_protions(pdb_file_name='1bbd',chain_ID='H',limit=113)
  pdb_hierarchy_var_L,pdb_hierarchy_const_L = get_pdb_protions(pdb_file_name='1bbd',chain_ID='L',limit=103)
  sites_var_heavy,sites_const_heavy,selection_var_H,selection_const_H,pdb_hierarchy = get_sites_protions(pdb_file_name='1bbd',chain_ID='H',limit=113)
  sites_var_light,sites_const_light,selection_var_L,selection_const_L,pdb_hierarchy = get_sites_protions(pdb_file_name='1bbd',chain_ID='L',limit=107)


def get_file(fn):
  if not os.path.isfile(fn + '.pdb'):
    sf_fn = fetch.get_pdb (fn,data_type='pdb',mirror='rcsb',log=sys.stdout)


def set_work_folder():
  if sys.platform.startswith('win'):
    wrokpath = r'C:\Phenix\Dev\Work\work\FAB elbow angle'
  else:
    workpath = '~youval/Work/work/FAB'
  os.chdir(wrokpath)

def check_pdb_file_name(file_name):
  tmp = os.path.basename(file_name)
  tmp = tmp.split('.')
  assert len(tmp[0])==4
  if len(tmp)==2 and tmp[1]=='pdb':
    fn = file_name
  else:
    fn = tmp[0] + '.pdb'
  return fn

def get_pdb_protions(pdb_file_name,chain_ID,limit):
  '''(str,str,int) -> pdb_hierarchy,pdb_hierarchy
  create a pdb_hierarchy from pdb_file_name, takes a chain,
  with chain_ID, and devides it to two, at the limit position.


  Example:
  >>>pdb_heavy_var, pdb_heavy_const = get_pdb_protions(pdb_file_name='1bbd',chain_ID='H',limit=113)

  Arguments:
  ----------
  pdb_file_name : pdb file name
  chain_ID : PDB chain ID
  limit : the cutoff between the Variable and the Constant FAB domains

  Returns:
  --------
  pdb_hierarchy_var : 0 to limit portion of the pdb_hierarchy
  pdb_hierarchy_const : limit to end portion of the pdb_hierarchy
  '''
  fn = check_pdb_file_name(pdb_file_name)
  pdb_hierarchy = pdb.hierarchy.input(file_name=fn).hierarchy
  chains = pdb_hierarchy.models()[0].chains()
  chain_names = [x.id for x in chains]
  # check if correct name is given and find the correct chain
  if chain_ID in chain_names:
    # get the chain residue group, to check its length
    indx = chain_names.index(chain_ID)
    # get the chain with the chain_ID and index indx
    chain_res = pdb_hierarchy.models()[0].chains()[indx].residue_groups()
    chain_length = len(chain_res)
    selection_var = 'chain {0} resseq {1}:{2}'.format(chain_ID,1,limit)
    selection_const = 'chain {0} resseq {1}:{2}'.format(chain_ID,limit,chain_length)
    pdb_hierarchy_var = pdb_hierarchy.atom_selection_cache().selection(selection_var)
    pdb_hierarchy_const = pdb_hierarchy.atom_selection_cache().selection(selection_const)
  else:
    raise Sorry('The is not chain with {0} name in {1}'.format(chain_ID,fn))

  return pdb_hierarchy_var,pdb_hierarchy_const

def get_sites_protions(pdb_file_name,chain_ID,limit):
  '''(str,str,int) -> pdb_hierarchy,pdb_hierarchy
  create a pdb_hierarchy from pdb_file_name, takes a chain,
  with chain_ID, and devides it to two, at the limit position.


  Example:
  >>>pdb_heavy_var, pdb_heavy_const = get_pdb_protions(pdb_file_name='1bbd',chain_ID='H',limit=113)

  Arguments:
  ----------
  pdb_file_name : pdb file name
  chain_ID : PDB chain ID
  limit : the cutoff between the Variable and the Constant FAB domains

  Returns:
  --------
  sites_var,sites_const
  sites_var : atoms sites of 0 to limit portion of the pdb_hierarchy
  sites_const : atoms sites of limit to end portion of the pdb_hierarchy
  '''
  fn = check_pdb_file_name(pdb_file_name)
  pdb_inp = pdb.input(file_name=fn)
  pdb_hierarchy = pdb_inp.construct_hierarchy()
  xrs = pdb_inp.xray_structure_simple()
  sites_cart_all = xrs.sites_cart()
  chains = pdb_hierarchy.models()[0].chains()
  chain_names = [x.id for x in chains]
  # check if correct name is given and find the correct chain
  if chain_ID in chain_names:
    # get the chain residue group, to check its length
    indx = chain_names.index(chain_ID)
    # get the chain with the chain_ID and index indx
    chain_res = pdb_hierarchy.models()[0].chains()[indx].residue_groups()
    chain_length = len(chain_res)
    selection_var = 'chain {0} resseq {1}:{2}'.format(chain_ID,1,limit)
    selection_const = 'chain {0} resseq {1}:{2}'.format(chain_ID,limit,chain_length)
    pdb_var_selection = pdb_hierarchy.atom_selection_cache().selection(string=selection_var)
    pdb_const_selection = pdb_hierarchy.atom_selection_cache().selection(string=selection_const)
    # get atoms
    params = superpose_pdbs.master_params.extract()
    x = superpose_pdbs.manager(
      params,
      log='log.txt',
      pdb_hierarchy_moving=pdb_const_selection,
      pdb_hierarchy_fixed=pdb_var_selection)

    sites_var = sites_cart_all.select(pdb_var_selection)
    sites_const = sites_cart_all.select(pdb_const_selection)

    chain = chain_res = pdb_hierarchy.models()[0].chains()[indx]

    x1 = superpose_pdbs.manager.extract_sequence_and_sites_chain(chain,'resseq 0:113')
    x2 = superpose_pdbs.manager.extract_sequence_and_sites_chain(chain,'resseq 113:190')
  else:
    raise Sorry('The is not chain with {0} name in {1}'.format(chain_ID,fn))

  return sites_var,sites_const,selection_var,selection_const

class FAB_elbow_angle(object):
  def __init__(self,
               file_nam,
               chain_name_light,
               chain_name_heavy,
               limit_light=107,
               limit_heavy=113):
    '''Get elbow angle for Fragment antigen-binding (Fab)'''


if __name__=='__main__':
  fn = '1bbd'
  set_work_folder()
  get_file(fn)
  run(fn)
  print 'Done'
