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
  pdb_hierarchy = pdb.hierarchy.input(file_name=fn).hierarchy
  #select_str = 'chain H and resseq 1:end'
  select_str = 'chain L not HETATM'
  #select_str = 'chain H and resseq 0:102'
  pdb_selection_bool = pdb_hierarchy.atom_selection_cache().selection(select_str)
  selected_chain_portion  = pdb_hierarchy.select(pdb_selection_bool)
  selected_chain_portion.atoms()[0].xyz

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
