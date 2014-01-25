from __future__ import division
from FAB.FAB_elbow_angle import FAB_elbow_angle
from iotbx.pdb import fetch
import os,sys


def run(fn):
  ''' '''
  fab = FAB_elbow_angle(
    pdb_file_name=fn,
    chain_ID_light='L',
    chain_ID_heavy='H',
    limit_light=107,
    limit_heavy=113)

  print 'ok'


def get_file(fn):
  '''If file is not in local directory, it will me fetched from rcsb website'''
  if not os.path.isfile(fn + '.pdb'):
    pdb_fn = fetch.get_pdb (fn,data_type='pdb',mirror='rcsb',log=sys.stdout)



def set_work_folder():
  if sys.platform.startswith('win'):
    wrokpath = r'C:\Phenix\Dev\Work\work\FAB'
  else:
    workpath = '~youval/Work/work/FAB'
  os.chdir(wrokpath)




if __name__=='__main__':
  fn = '1bbd'
  set_work_folder()
  get_file(fn)
  run(fn)
  print 'Done'