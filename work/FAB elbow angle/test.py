from __future__ import division
from phenix.command_line import superpose_pdbs
from iotbx.pdb import fetch
from libtbx import easy_run
import os,sys


def get_file(fn):
  if not os.path.isfile(fn + '.pdb'):
    sf_fn = fetch.get_pdb (fn,data_type='pdb',mirror='rcsb',log=sys.stdout)


def set_work_folder():
  if sys.platform.startswith('win'):
    wrokpath = r'C:\Phenix\Dev\Work\work\FAB elbow angle'
  else:
    workpath = '~youval/Work/work/FAB'
  os.chdir(wrokpath)

if __name__=='__main__':
  set_work_folder()
  get_file('1bbd')
  # Test if FAB

  # Devide to



  print 'Done'
