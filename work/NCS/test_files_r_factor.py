from __future__ import division
from  misc_scripts.r_factor_calc import *
from libtbx import easy_run
import os

def run(pdb_file,cif_file):
  '''
  Test of a single file
  '''
  t = r_factor_calc(
    [pdb_file,cif_file],
    eps=1e-2,
    strOut=True,fromRCSB=True)
  print t

if __name__=='__main__':
  os.chdir(r'C:\Phenix\Dev\Work\work\NCS\junk\pdb_test')
  pdb_id = '1b35'
  run(pdb_file=pdb_id+'.pdb', cif_file=pdb_id+'-sf.cif')


