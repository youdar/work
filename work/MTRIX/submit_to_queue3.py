from __future__ import division
from libtbx.command_line import easy_qsub
from misc_scripts import helpers
from  misc_scripts.r_factor_calc import *
from iotbx import pdb
import cPickle as pickle
import os

'''
Submit to queue all files with good MTRIX records that also havestracture factor files
'''

def run():
  # set working environment
  #phenix_source = "/net/chevy/raid1/youval/Work_chevy/phenix_build/setpaths.csh"
  #where_to_run_dir = "/net/cci-filer2/raid1/home/youval/Work/work/junk/queue_job_3"
  phenix_source = r"c:\Phenix\Dev\phenix_build\setpaths.csh"
  #os.chdir(where_to_run_dir)
  # collect files to work on
  files = get_file_names()
  #
  commands = []
  # Set command path
  com_path = '/net/chevy/raid1/youval/Work_chevy/phenix_sources/misc_scripts/r_factor_calc.py'
  #
  for [pdb_file,sf_file,file_name] in files:
    #print file_name
    com_options = '--strOut=True --eps=0.01 --fromRCSB=False'.format(file_name)
    outString = '{0} {1} {2} {3} &> log_{4}'.format(com_path,pdb_file,sf_file,com_options,file_name)
    commands.append(
      "python {}".format(outString))

  easy_qsub.run(
    phenix_source = phenix_source,
    where         = where_to_run_dir,
    qsub_cmd      = 'qsub -q all.q@morse',
    commands      = commands,
    size_of_chunks= 300) # because jobs are quick

def get_file_names():
  '''
  returns the list of files that need to be processed
  '''
  #data_dir = '/net/cci-filer2/raid1/home/youval/Work/work/MTRIX/Data'
  data_dir = r'c:\Phenix\Dev\Work\work\MTRIX\Data'
  results = []
  # all 157 files with good MTRIX also have structure factors files
  # good MTRIX is tested with eps = 0.01
  files_with_good_MTRIX = set(pickle.load(open(os.path.join(data_dir,'files_with_good_MTRIX'),'r')))
  # file location dictionary
  good_MTRIX_pdb_files = pickle.load(open(os.path.join(data_dir,'dict_good_MTRIX_pdb_files'),'r'))
  structure_factors_files = pickle.load(open(os.path.join(data_dir,'dict_structure_factors_files'),'r'))
  print 'File names and location dictionaries are loaded...'
  # Just to explore data
  mtrix_not_included_in_pdb = set(open(os.path.join(data_dir,'mtrix_not_included_in_pdb.txt')).read().splitlines())

  # look at files that were processed
  files_with_problems = pickle.load(open(os.path.join(data_dir,'files_with_problems'),"r"))
  files_with_problems = {x[0] for x in files_with_problems}

  Collect_tested_files = pickle.load(open(os.path.join(data_dir,'Collect_tested_files'),"r"))
  Collect_tested_files = {x[0] for x in Collect_tested_files}

  file_set = (files_with_good_MTRIX - Collect_tested_files)
  print 'Number of files from the 157 that have issues: {}'.format(len(mtrix_not_included_in_pdb - Collect_tested_files))
  print 'number of files_with_good_MTRIX: {}'.format(len(files_with_good_MTRIX))
  print 'number of file records in Collect_tested_files: {}'.format(len(Collect_tested_files))
  print 'number of files to run: {}'.format(len(file_set))

  #file_set = set(['3dpr', '1o4z', '1bxb', '2bht', '4a5p', '4a5m', '2x9m',
                  #'4bl4', '2x9i', '4gme', '3lob', '2vvq', '1g31', '2wy3',
                  #'4iq1', '2wy6', '2w0c', '2uwa'])


  for file_name in file_set:
      pdb_file = good_MTRIX_pdb_files[file_name]
      sf_file = structure_factors_files[file_name]
      results.append([pdb_file,sf_file,file_name])
  return results

if (__name__ == "__main__"):
  run()
