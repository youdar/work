from __future__ import division
from libtbx.command_line import easy_qsub
from misc_scripts import helpers
from  misc_scripts.r_factor_calc import *
from iotbx import pdb
import cPickle as pickle
import random
import os,sys


def run():
  # convert the path to python format
  directory_path = '/net/cci/youval/Work/work/Clashes/Data'
  directory_path = os.path.realpath(directory_path)
  os.chdir('/net/cci/youval/Work/work/Clashes/wtest')
  # pdb_clash_score_and_name = list([score_with_hydrogen,score_without_hydrogen,experment_type,file_name]...)
  pdb_clash_score_and_name = pickle.load(open(os.path.join(directory_path,'pdb_clash_score_and_name'),'r'))
  pdb_clash_score_and_name.sort()
  # only look at files with clashscore at a certain range
  file_list = []
  while 1:
    x = pdb_clash_score_and_name.pop(0)
    if x[0] > 50: break
    file_list.append(x[3])
  # selecet n files randomly
  #n = 1
  #files = random.sample(file_list,n)
  files = file_list
  # collect files full path, path in the PDB_MIRROR_PDB, to work on
  #files_list = get_file_names(files)
  files_list = [[x+'.pdb',x] for x in files]
  # set working environment
  phenix_source = "/net/chevy/raid1/youval/Work_chevy/phenix_build/setpaths.csh"
  where_to_run_dir = "/net/cci/youval/Work/work/Clashes/queue_clash_compare_12_11_2013"
  os.chdir(where_to_run_dir)
  commands = []
  cntr = 0
  # Set command path
  com_path = '/net/cci/youval/Work/work/Clashes/Test_internal_clashscore.py'
  #
  for [pdb_file,file_name] in files_list:
    outString = '{0} {1} &> log_{2}'.format(com_path,pdb_file,file_name)
    commands.append(
      "python {}".format(outString))


  easy_qsub.run(
    phenix_source = phenix_source,
    where         = where_to_run_dir,
    commands      = commands,
    #qsub_cmd      = 'qsub -q all.q@rebus',
    #qsub_cmd      = 'qsub -q all.q@rebus',
    qsub_cmd      = 'qsub -q all.q@morse',
    size_of_chunks= 300) # because jobs are quick

def get_file_names(file_list):
  '''
  returns the list of files that need to be processed
  '''
  #os.chdir('/net/cci-filer2/raid1/home/youval/Work/work')
  results = []
  for f in file_list:
    file_path = find_full_path(f)
    results.append([file_path,f])
  return results

def find_full_path(file_name):
  '''(str) -> str
  returns the full path a pdb file name in  LBL PDB_MIRROR_PDB
  '''
  # change file name to the form in the PDB_MIRROR_PDB
  f = 'pdb{}.ent.gz'.format(file_name)
  full_path = []
  pdb_dir = os.environ["PDB_MIRROR_PDB"]
  for root, _, files in os.walk(pdb_dir):
    if f in files:
      full_path = os.path.join(root, f)
      break
  return full_path


if (__name__ == "__main__"):
  run()