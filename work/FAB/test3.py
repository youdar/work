from __future__ import division
from phenix.command_line import superpose_pdbs
from iotbx import pdb
from iotbx.pdb import fetch
from libtbx.utils import Sorry
from libtbx.utils import null_out
import shutil
import tempfile
import os,sys

def run(fn):
  ''' '''
  pdb_hierarchy = pdb.hierarchy.input(file_name=fn).hierarchy
  ss = 'pepnames and (name ca or name n or name c) and altloc " "'
  select_str = 'chain L and resseq 104:end and %s' % ss
  select_str = 'chain H and resseq 114:end and pepnames and (name ca or name n or name c) and altloc " "'
  pdb_selection_bool = pdb_hierarchy.atom_selection_cache().selection(select_str)
  const_L  = pdb_hierarchy.select(pdb_selection_bool)


  ss = 'pepnames and (name ca or name n or name c) and altloc " "'
  select_str = 'chain H and resseq 114:end and %s' % ss
  pdb_selection_bool = pdb_hierarchy.atom_selection_cache().selection(select_str)
  const_H  = pdb_hierarchy.select(pdb_selection_bool)

  params = superpose_pdbs.master_params.extract()



  x = superpose_pdbs.manager(
    params,
    log=null_out(),
    write_output=False,
    save_lsq_fit_obj=True,
    pdb_hierarchy_fixed=const_H,
    pdb_hierarchy_moving=const_L)


  print 'ok'



def get_file(fn):
  if len(fn)==4: fn = fn + '.pdb'
  if not os.path.isfile(fn):
    fn = fetch.get_pdb (fn,data_type='pdb',mirror='rcsb',log=sys.stdout)
  return fn

def set_work_folder():
  if sys.platform.startswith('win'):
    wrokpath = r'C:\Phenix\Dev\Work\work\FAB'
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

if __name__=='__main__':
  fn = '1bbd'
  set_work_folder()
  fn = get_file(fn)
  run(fn)
  print 'Done'
