from __future__ import division
from  misc_scripts.r_factor_calc import *
from libtbx import easy_run
import cPickle as pickle
import os

def run(file_name):
  '''
  Test of a single wile
  '''
  good_MTRIX_pdb_files = pickle.load(open('dict_good_MTRIX_pdb_files','r'))
  good_BIOMT_pdb_files = pickle.load(open('dict_good_BIOMT_pdb_files','r'))
  structure_factors_files = pickle.load(open('dict_structure_factors_files','r'))
  MTRIX_with_Straucture_Factor = pickle.load(open('MTRIX_with_Straucture_Factor_file_list','r'))

  print 'working on {}'.format(file_name)
  pdb_file = good_MTRIX_pdb_files[file_name]
  sf_file = structure_factors_files[file_name]
  # calculate the precent of difference of R-work reconstructed vs mtz data
  #t = r_factor_calc([pdb_file,sf_file],eps=1e-3)
  #print t
  t = r_factor_calc([pdb_file,sf_file],eps=1e-3,file_name=file_name,strOut=True)
  print t

def get_file(fn):
  os.chdir('/net/cci-filer2/raid1/home/youval/Work/work/junk')
  cmd = "phenix.fetch_pdb {}".format(fn)
  r = easy_run.go(cmd)
  os.chdir('/net/cci-filer2/raid1/home/youval/Work/work')


if __name__=='__main__':
  # move to working directory
  os.chdir('/net/cci-filer2/raid1/home/youval/Work/work')
  #file_name = '4kn2' # have both IOBS and FOBS
  #file_name = '4aun'  # have issues running phenix.cif_as_mtz
  #file_name = '2wws'
  run('3nby')