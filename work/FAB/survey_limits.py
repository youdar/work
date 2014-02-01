from __future__ import division
from mmtbx.utils.FAB_elbow_angle import FAB_elbow_angle
from libtbx.utils import null_out
from iotbx.pdb import fetch
import itertools
from pprint import pprint
import os,sys

'''
Exploring how selection of limit between variable and constant chains
effect the Fragment antigen-binding (Fab) angle

There are two limits, one for the heavy chain and one for the light one
I will sample +- 3 rediues
'''


class explor_fab(object):
  '''
  '''
  def __init__(self,file_name):
    # set test path
    osType = sys.platform
    if osType.startswith('win'):
        os.chdir(r'C:\Phenix\Dev\Work\work\FAB\junk')
    else:
        os.chdir('/net/cci-filer2/raid1/home/youval/Work/work/FAB/junk')
    if not os.path.isfile(file_name):
      fn = fetch.get_pdb ('1bbd',data_type='pdb',mirror='rcsb',log=null_out())
    self.limits_delta = []
    self.file_name = fn
    tmp = os.path.basename(file_name)
    self.pdb_code = tmp.split('.')[0]

  def set_test_range(self,n):
    self.limits_delta = [x for x in itertools.permutations(range(-n,n+1),2)]
    self.limits_delta.append((0,0))

  def get_FAB_angle(self,limit_light=107,limit_heavy=113):
    fab = FAB_elbow_angle(pdb_file_name=self.file_name,limit_light=limit_light,limit_heavy=limit_heavy)
    return fab.FAB_elbow_angle

  def calc_fab(self,range=0,limit_light=107,limit_heavy=113):
    ''''''
    self.set_test_range(range)
    self.fab_angles = []
    for (l_L,l_H) in self.limits_delta:
      angel = self.get_FAB_angle(
        limit_light=limit_light + l_L,
        limit_heavy=limit_heavy + l_H)
      self.fab_angles.append([limit_light + l_L, limit_heavy + l_H, int(angel)])

  def fab_angle_print(self):
    ''''''
    print 'PDB code {0}, limits survey results'.format(self.pdb_code)
    print '------------------------------------------'
    pprint(self.fab_angles)



if __name__=='__main__':

  for name in ['1bbd','7fab','1dba','1plg','1nl0']:
    t = explor_fab(file_name='1bbd')
    t.calc_fab(range=2)
    t.fab_angle_print()
    t.fab_angle_print()
