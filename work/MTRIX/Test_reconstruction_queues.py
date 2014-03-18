from __future__ import division
from  misc_scripts.r_factor_calc import *
from  iotbx.pdb.multimer_reconstruction import multimer
from multiprocessing import Process, Queue, Lock
from iotbx import pdb
import cPickle as pickle
import os


'''
Read list of pdb files names with more than one good BIOMT records
Read list of pdb files names with more than one good MTRIX records

Get coresponding structure factor files

@author: Youval Dar
'''

def Call_function(queue,lock):
  '''
  Collect the results from the parallel process and write them into files

  Collect_tested_files : the file that lists all pdb files that were checked
  files_with_problems : the files with all pdb files that had processing issues and r>=1

  results are stored in files in the following format:

  file_name1:r1:msg1
  file_name2:r2:msg2
  file_name3:r3:msg3
  .
  .

  the r is the result of r_factor_calc.py
  msg is the error or problem with the test, if there is one

  '''
  # append results from this run
  f = open('/net/cci-filer2/raid1/home/youval/Work/work/Collect_tested_files','a')
  g = open('/net/cci-filer2/raid1/home/youval/Work/work/files_with_problems','a')

  while True:
    # get items from the queue
    x = queue.get()
    # check if queue is empty
    if x == 'DONE':
      f.close()
      g.close()
      break
    # if we are not DONE
    # calculate the precent of difference of R-work reconstructed vs mtz data
    [pdb_file,sf_file,file_name] = x
    print 'Processing file {}'.format(file_name)
    try:
      r = r_factor_calc([pdb_file,sf_file],eps=2e-3,fromRCSB=False)
      msg = 'OK'
    except Sorry as e:
      r = 100
      msg = e.message
    except TypeError as e:
      r = 100
      msg = e.message
    # Write results to file
    outString = '{0}:{1}:{2}\n'.format(file_name,r,msg)
    lock.acquire()
    if r<1:
      f.write(outString)
      print outString
    else:
      g.write(outString)
      print outString
    lock.release()




def make_dict(index_file_name,data_dir=''):
  '''
  Read all file names from PBP mirror folder, structure factor files
  or other file containing file names

  for PDB fils check the correct folder using os.environ["PDB_MIRROR_PDB"]
  and the file name is 'INDEX'
  for structure factor files use os.environ["PDB_MIRROR_STRUCTURE_FACTORS"]
  and the file name 'INDEX'

  input:
  data_dir : the directory containing a file with the names of files we want to extract
  index_file_name : file names list

  Output:
  a dictionary
  {file_name: file_path,...}
  '''
  file_names = open(os.path.join(data_dir, index_file_name), "r").readlines()
  result = {}
  for file_path in file_names:
    # file_path looks like: '/net/chevy/raid1/pdb_mirror/pdb/00/pdb100d.ent.gz'
    # file_name should look like 'pdb100d'
    file_path = file_path.strip()
    file_name = file_path.split('/')[-1]
    file_name = file_name.split('.')[0]
    if file_name.startswith('pdb'):
      # pdb file names in INDEX are like 'pdb2vpf', the file_name should be '2vpf'
      file_name = file_name[3:]
    elif file_name.startswith('r'):
      # structure factor file names in INDEX are like 'r1oqjsf', it should be '1oqj'
      file_name = file_name[1:-2]
    else:
      print 'File namming problems!!!'
      print file_name
      break
    if file_path.startswith('/net'):
      result[file_name] = file_path
    else:
      result[file_name] = data_dir+file_path
  return result


def run():
  '''
  good_MTRIX_pdb_files, good_BIOMT_pdb_files and structure_factors_files
  are dictionaries. the keys are pdb record name and the values are the
  appropriate file full path
  '''

  # If you already have the dictionaries use:
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
  # Test of all files in MTRIX_with_Straucture_Factor
  # collect all good results and save them on a file so that
  # not to repeat them
  #tested_files = open('Collect_tested_files',"r").readlines()
  tested_files = open('Collect_tested_files',"r").read().splitlines()
  files_with_problems = open('files_with_problems',"r").read().splitlines()
  # Clean the remarks - use only protein name
  files_with_problems = [x[:4] for x in files_with_problems]
  tested_files = [x[:4] for x in tested_files]
  # start queue
  queue = Queue()
  lock = Lock()
  reader_p = Process(target=Call_function, args=(queue,lock))
  reader_p.daemon = True
  reader_p.start()
  for file_name in MTRIX_with_Straucture_Factor:
    if (file_name not in tested_files) and (file_name not in files_with_problems):
      pdb_file = good_MTRIX_pdb_files[file_name]
      sf_file = structure_factors_files[file_name]
      queue.put([pdb_file,sf_file,file_name])
  queue.put('DONE')
  # close the parallel process
  reader_p.join()		# wait for the reader to finish

  # Analyze the results and add them to dictionaries
  tested_files = open('Collect_tested_files',"r").read().splitlines()
  for x in tested_files:
    [file_name,r,msg] = x.split(':')
    r = float(r)
    reconstruction_test_dict[file_name] = r
    reconstruction_test_list.append(r)
  # save the results
  pickle.dump(reconstruction_test_dict,open('reconstruction_test_dict','w'))
  pickle.dump(reconstruction_test_list,open('reconstruction_test_list','w'))

  print 'Done...'


if __name__=='__main__':
  # move to working directory
  os.chdir('/net/cci-filer2/raid1/home/youval/Work/work')
  #os.chdir('c:\\Phenix\\Dev\\Work\\work')
  # check how many processors are available
  run()
