from __future__ import division
from libtbx.command_line import easy_qsub
from misc_scripts import helpers
from  misc_scripts.r_factor_calc import *
from iotbx import pdb
import cPickle as pickle
import os

def get_file_names():
  '''
  returns the list of files that need to be processed
  '''
  os.chdir('/net/cci-filer2/raid1/home/youval/Work/work')
  results = []
  good_MTRIX_pdb_files = pickle.load(open('dict_good_MTRIX_pdb_files','r'))
  #good_BIOMT_pdb_files = pickle.load(open('dict_good_BIOMT_pdb_files','r'))
  structure_factors_files = pickle.load(open('dict_structure_factors_files','r'))
  MTRIX_with_Straucture_Factor = pickle.load(open('MTRIX_with_Straucture_Factor_file_list','r'))
  # Load previous results
  #reconstruction_test_dict = pickle.load(open('reconstruction_test_dict','r'))
  #reconstruction_test_list = pickle.load(open('reconstruction_test_list','r'))
  print 'Dictionaries are loaded...'
  # iterate over file and calculate qulity of R-work of reconstructed pdb file
  tested_files = open('Collect_tested_files',"r").read().splitlines()
  files_with_problems = open('files_with_problems',"r").read().splitlines()
  files_with_problems = [x[:4] for x in files_with_problems]
  tested_files = [x[:4] for x in tested_files]
  print 'Previous processed file list loaded'
  for file_name in MTRIX_with_Straucture_Factor:
      if (file_name not in tested_files) and (file_name not in files_with_problems):
        pdb_file = good_MTRIX_pdb_files[file_name]
        sf_file = structure_factors_files[file_name]
        results.append([pdb_file,sf_file,file_name])
  return results

def run():
  # set working environment
  os.chdir('/net/cci/youval/Work/work/MTRIX')
  phenix_source = "/net/chevy/raid1/youval/Work_chevy/phenix_build/setpaths.csh"
  where_to_run_dir = "/net/cci-filer2/raid1/home/youval/Work/work/queue_job"
  # collect files to work on
  files = get_file_names()
  #
  commands = []
  #cntr = 0
  # Set command path
  com_path = '/net/chevy/raid1/youval/Work_chevy/phenix_sources/misc_scripts/r_factor_calc.py'

  for [pdb_file,sf_file,file_name] in files:
    #print file_name
    com_options = '--strOut=True'.format(file_name)
    outString = '{0} {1} {2} &> log_{1}'.format(com_path,file_name,com_options)
    commands.append(
      "python {}".format(outString))
    #cntr += 1
    #if(cntr == 50000): break
  easy_qsub.run(
    phenix_source = phenix_source,
    where         = where_to_run_dir,
    commands      = commands,
    size_of_chunks= 200) # because jobs are quick

if (__name__ == "__main__"):
  run()
