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

# global variable for parallel process results collection
results = []

def collect_results(x):
  '''
  Collect the results from the parallel process
  '''
  results.append(x)

def full_pdb_file_paths(file_type=1):
  '''
  read all names of
  1: PDB files
  2: structure factor files

  return:
  {file_name: file_path,...}
  '''
  if file_type == 1:	# PDB files
    files_dir = os.environ["PDB_MIRROR_PDB"]
  else:			# structure factor files
    files_dir = '/net/cci/pdb_mirror/structure_factors'

  file_names = open(os.path.join(files_dir, "INDEX"), "r").readlines()
  result = {}
  for file_path in file_names:
    file_path = os.path.join(files_dir, file_path.strip())
    file_name = file_path.split('/')[-1]
    file_name = file_name.split('.')[0]
    result[file_name] = file_path
  return result

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

def start_multiprocessing():
  while True:
    try:
      np = int(raw_input('Number of processors available is {}. How many would you like to use? '.format(mp.cpu_count())))
      if np>mp.cpu_count():
        raise Exception
      break
    except ValueError:
      print 'Please enter an integer \n'
    except Exception:
      print 'Please choose smaller number of processors \n'
  # set number of CPUs
  p = mp.Pool(processes=np)
  return p
def run(recon_test=False,build_new_dictinaries=False):
  '''
  good_MTRIX_pdb_files, good_BIOMT_pdb_files and structure_factors_files
  are dictionaries. the keys are pdb record name and the values are the
  appropriate file full path
  '''

  if build_new_dictinaries:
    # Do the following only if there were changes to the files lists
    good_MTRIX_pdb_files = make_dict('mtrix_ok_run.txt')
    good_BIOMT_pdb_files = make_dict('biomt_ok_run.txt')
    structure_factors_files = make_dict('INDEX','/net/cci/pdb_mirror/structure_factors/')
    pickle.dump(good_MTRIX_pdb_files,open('dict_good_MTRIX_pdb_files','w'))
    pickle.dump(good_BIOMT_pdb_files,open('dict_good_BIOMT_pdb_files','w'))
    pickle.dump(structure_factors_files,open('dict_structure_factors_files','w'))
    print 'Dictionaries Created...'
  else:
    # If you already have the dictionaries use:
    good_MTRIX_pdb_files = pickle.load(open('dict_good_MTRIX_pdb_files','r'))
    good_BIOMT_pdb_files = pickle.load(open('dict_good_BIOMT_pdb_files','r'))
    structure_factors_files = pickle.load(open('dict_structure_factors_files','r'))
    MTRIX_with_Straucture_Factor = pickle.load(open('MTRIX_with_Straucture_Factor_file_list','r'))
    print 'Dictionaries are loaded...'

  # When changing the file lists
  if build_new_dictinaries:
    # When changing the file lists
    MTRIX_with_Straucture_Factor = []
    for x in good_MTRIX_pdb_files:
      if structure_factors_files.has_key(x):
        MTRIX_with_Straucture_Factor.extend([x,good_MTRIX_pdb_files[x],structure_factors_files[x]])

    l = len(good_MTRIX_pdb_files)
    i = len(MTRIX_with_Straucture_Factor)
    print 'The number of both structure factors and good MTRIX with the same name: {} from a total of {}'.format(i,l)
    pickle.dump(MTRIX_with_Straucture_Factor,open('MTRIX_with_Straucture_Factor_file_list','w'))

    i = 0
    for x in good_BIOMT_pdb_files:
      if structure_factors_files.has_key(x):
        i += 1
    l = len(good_BIOMT_pdb_files)
    print 'The number of both structure factors and good BIOMT with the same name: {} from a total of {}'.format(i,l)

    #f1 = open('dict_good_MTRIX_pdb_files.txt','w')
    #f1.writelines([x+'\n' for x in dict_good_MTRIX_pdb_files])
    #f1.close()
  # run test - compare r-work fromreconstructed pdb file to that of the mtz data
  if recon_test:
    p = start_multiprocessing()
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
    # append results from this run
    f = open('Collect_tested_files',"a")
    g = open('files_with_problems',"a")
    for file_name in MTRIX_with_Straucture_Factor:
      if (file_name not in tested_files) and (file_name not in files_with_problems):
        print file_name
        pdb_file = good_MTRIX_pdb_files[file_name]
        sf_file = structure_factors_files[file_name]
        # calculate the precent of difference of R-work reconstructed vs mtz data
        p.apply_async(r_factor_calc,([pdb_file,sf_file],),
                      {'eps':2e-3,'file_name':file_name,'strOut':True},
                      callback=collect_results)
    # close the parallel process
    p.close()
    p.join()
    # The multiprocess does not catch all the raised errors and sorrys
    # need to run the regular version after this is done
    print '*'*50
    print 'Done with calculating'
    print '*'*50
    for x in results:
      # x is of the form 'pdbname:r_score'
      [file_name,r] = x.split(':')
      r = float(r)
      reconstruction_test_dict[file_name] = r
      reconstruction_test_list.append(r)
      if r<1:
        f.write(x + '\n')
      else:
        g.write(x + '\n')

    f.close()
    g.close()


    ## Test of a single wile
    ##file_name = '4kn2' # have both IOBS and FOBS
    #file_name = '4aun'  # have issues running phenix.cif_as_mtz
    #file_name = '2wws'
    #print file_name
    #pdb_file = good_MTRIX_pdb_files[file_name]
    #sf_file = structure_factors_files[file_name]
    ## calculate the precent of difference of R-work reconstructed vs mtz data
    #t = r_factor_calc([pdb_file,sf_file],eps=1e-3)
    #reconstruction_test_dict[file_name] = t
    #reconstruction_test_list.append(t)

    # save the results
    pickle.dump(reconstruction_test_dict,open('reconstruction_test_dict','w'))
    pickle.dump(reconstruction_test_list,open('reconstruction_test_list','w'))

  print 'Done...'


if __name__=='__main__':
  # move to working directory
  os.chdir('/net/cci-filer2/raid1/home/youval/Work/work')
  #os.chdir('c:\\Phenix\\Dev\\Work\\work')
  # check how many processors are available
  run(recon_test=True)
