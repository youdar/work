from __future__ import division
import os
from misc_scripts.strict_ncs import test_ncs_refinement as tst

__author__ = 'Youval'


def run(fn):
  tst.run([fn])
  print 'Done...'


if __name__ == '__main__':
  os.chdir(r'C:\Phenix\Dev\Work\work\NCS\junk')
  fn = '4m40.pdb'
  run(fn)
