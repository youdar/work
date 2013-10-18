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
  pdb_dir = os.environ["PDB_MIRROR_PDB"]
  pdb_files = open(os.path.join(pdb_dir, "INDEX"), "r").read().splitlines()
  for f in pdb_files:
    file_name = f.split('/')[-1]
    file_name = file_name.split('.')[0][3:7]
    results.append([os.path.join(pdb_dir,f),file_name])
  return results

def run():
  # set working environment
  phenix_source = "/net/chevy/raid1/youval/Work_chevy/phenix_build/setpaths.csh"
  where_to_run_dir = "/net/cci/youval/Work/work/Clashes/queue_clash"
  # collect files to work on
  files = get_file_names()
  #
  commands = []
  cntr = 0
  # Set command path
  com_path = '/net/cci/youval/Work/work/Clashes/test_clashes.py'
  #
  for [pdb_file,file_name] in files:
    outString = '{0} {1} &> log_{2}'.format(com_path,pdb_file,file_name)
    commands.append(
      "python {}".format(outString))


  easy_qsub.run(
    phenix_source = phenix_source,
    where         = where_to_run_dir,
    commands      = commands,
    size_of_chunks= 200) # because jobs are quick

if (__name__ == "__main__"):
  run()
