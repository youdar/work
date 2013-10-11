from __future__ import division
from  misc_scripts.r_factor_calc import *
from  iotbx.pdb.multimer_reconstruction import multimer
import multiprocessing as mp
from iotbx import pdb
import cPickle as pickle
import os


'''
Read list of pdb files names with more than one good BIOMT records
Read list of pdb files names with more than one good MTRIX records

Get coresponding structure factor files

@author: Youval Dar
'''

def run():
  '''
  good_MTRIX_pdb_files, good_BIOMT_pdb_files and structure_factors_files
  are dictionaries. the keys are pdb record name and the values are the
  appropriate file full path
  '''

  good_MTRIX_pdb_files = pickle.load(open('dict_good_MTRIX_pdb_files','r'))
  good_BIOMT_pdb_files = pickle.load(open('dict_good_BIOMT_pdb_files','r'))
  structure_factors_files = pickle.load(open('dict_structure_factors_files','r'))
  MTRIX_with_Straucture_Factor = pickle.load(open('MTRIX_with_Straucture_Factor_file_list','r'))
  print 'Dictionaries are loaded...'
  # run test - compare r-work fromreconstructed pdb file to that of the mtz data
  print '*'*50
  print 'Start testing MTRIX reconstruction testing'
  print '*'*50
  # Load previous results
  reconstruction_test_dict = pickle.load(open('reconstruction_test_dict','r'))
  reconstruction_test_list = pickle.load(open('reconstruction_test_list','r'))
  # iterate over file and calculate qulity of R-work of reconstructed pdb file
  tested_files = open('Collect_tested_files',"r").read().splitlines()
  files_with_problems = open('files_with_problems',"r").read().splitlines()
  files_with_problems = [x[:4] for x in files_with_problems]
  f = open('Collect_tested_files',"a")
  g = open('files_with_problems',"a")
  # iterate over file and calculate qulity of R-work of reconstructed pdb file
  for file_name in MTRIX_with_Straucture_Factor:
    if (file_name not in tested_files) and (file_name not in files_with_problems):
      print file_name
      pdb_file = good_MTRIX_pdb_files[file_name]
      sf_file = structure_factors_files[file_name]
      # calculate the precent of difference of R-work reconstructed vs mtz data
      try:
        t = r_factor_calc([pdb_file,sf_file],eps=2e-3)
        msg = 'OK'
      except Sorry as e:
        msg = e.message
        t = 100
      except TypeError as e:
        msg = e.message
        t = 100
        outString = '{0}:{1}:{2}\n'.format(file_name,t,msg)
        if t<1:
          f.write(outString)
        else:
          g.write(outString)

  f.close()
  g.close()

  print 'Done...'


if __name__=='__main__':
  # move to working directory
  os.chdir('/net/cci-filer2/raid1/home/youval/Work/work')
  #os.chdir('c:\\Phenix\\Dev\\Work\\work')
  # check how many processors are available
  run()
