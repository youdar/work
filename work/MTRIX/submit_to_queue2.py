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
  os.chdir('/net/cci-filer2/raid1/home/youval/Work/work/junk')
  results = []
  good_MTRIX_pdb_files = pickle.load(open('dict_good_MTRIX_pdb_files','r'))
  #good_BIOMT_pdb_files = pickle.load(open('dict_good_BIOMT_pdb_files','r'))
  MTRIX_with_Straucture_Factor = pickle.load(open('MTRIX_with_Straucture_Factor_file_list','r'))
  files_with_mtrix_not_in_pdb = open('mtrix_not_included_in_pdb.txt',"r").read().splitlines()
  print 'Dictionaries are loaded...'
  # iterate over file and calculate qulity of R-work of reconstructed pdb file
  tested_files = open('Collect_tested_files',"r").read().splitlines()
  files_with_problems = open('files_with_problems',"r").read().splitlines()
  files_with_problems = [x[:4] for x in files_with_problems]
  tested_files = [x[:4] for x in tested_files]
  print 'Previous processed file list loaded'
  for file_name in MTRIX_with_Straucture_Factor:
    results.append([good_MTRIX_pdb_files[file_name],file_name])
  return results

def run():
  # set working environment
  phenix_source = "/net/chevy/raid1/youval/Work_chevy/phenix_build/setpaths.csh"
  where_to_run_dir = "/net/cci-filer2/raid1/home/youval/Work/work/queue_job_2"
  # collect files to work on
  files = get_file_names()
  #
  commands = []
  cntr = 0
  # Set command path
  com_path = '/net/cci-filer2/raid1/home/youval/Work/work/MTRIX/rotation_issues_list.py'
  #
  for [pdb_file,file_name] in files:
    outString = '{0} {1} &> log_{2}'.format(com_path,pdb_file,file_name)
    commands.append(
      "python {}".format(outString))


  easy_qsub.run(
    phenix_source = phenix_source,
    where         = where_to_run_dir,
    commands      = commands,
    size_of_chunks= 300) # because jobs are quick

if (__name__ == "__main__"):
  run()
