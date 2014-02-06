from __future__ import division
from mmtbx.utils.fab_elbow_angle import fab_elbow_angle
#from FAB.fab_elbow_angle_cross import fab_elbow_angle
from libtbx.utils import null_out
from iotbx.pdb import fetch
from iotbx import pdb
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
  def __init__(self,file_name,limit_L=107,limit_H=113):
    # set test path
    osType = sys.platform
    if osType.startswith('win'):
        os.chdir(r'C:\Phenix\Dev\Work\work\FAB\junk')
    else:
        os.chdir('/net/cci-filer2/raid1/home/youval/Work/work/FAB/junk')
    if not os.path.isfile(file_name):
      fn = fetch.get_pdb (file_name,data_type='pdb',mirror='rcsb',log=null_out())
    self.limits_delta = []
    self.file_name = fn
    self.limit_L = limit_L
    self.limit_H = limit_H
    tmp = os.path.basename(file_name)
    self.pdb_code = tmp.split('.')[0]

  def set_test_range(self,n):
    self.limits_delta = [x for x in itertools.permutations(range(-n,n+1),2)]
    self.limits_delta.append((0,0))

  def get_fab_angle(self,limit_light=None,limit_heavy=None):
    if not limit_light: limit_light=self.limit_L
    if not limit_heavy: limit_heavy=self.limit_H
    ph = pdb.input(file_name=self.file_name).construct_hierarchy()
    fab = fab_elbow_angle(pdb_hierarchy=ph,limit_light=limit_light,limit_heavy=limit_heavy)
    return fab.fab_elbow_angle

  def calc_fab(self,range=0,limit_light=107,limit_heavy=113):
    ''''''
    self.set_test_range(range)
    self.fab_angles = []
    for (l_L,l_H) in self.limits_delta:
      angel = self.get_fab_angle(
        limit_light=self.limit_L + l_L,
        limit_heavy=self.limit_H + l_H)
      self.fab_angles.append([self.limit_L + l_L, self.limit_H + l_H, int(angel)])

  def fab_angle_print(self):
    ''''''
    print 'PDB code {0}, limits survey results'.format(self.pdb_code)
    print 'Limit light, Limit heavy, Fab-angle'
    print '------------------------------------------'
    pprint(self.fab_angles)

  def write_to_file(self,file_name='survey_out.txt'):
    '''
    '''
    f = open(file_name,'w')
    f.write('PDB code {0}, limits survey results'.format(self.pdb_code))
    f.write('-'*30)
    for rec in self.fab_angles:
      f.write(rec)
    f.close()


def run(fab_list):
  for name,limit_L,limit_H in fab_list:
    t = explor_fab(file_name=name,limit_L=limit_L,limit_H=limit_H)
    t.calc_fab(range=2)
    t.fab_angle_print()

if __name__=='__main__':
  fab_list = sys.argv[1:]
  if not fab_list:
    fab_list = [('1bbd',114,118),('7fab',104,117),('1dba',107,113),('1plg',112,117),('1nl0',107,113)]
  run(fab_list)

