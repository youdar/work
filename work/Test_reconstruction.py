import os
from iotbx import pdb
import cPickle as pickle

'''
Read list of pdb files names with more than one good BIOMT records
Read list of pdb files names with more than one good MTRIX records

Get coresponding structure factor files

@author: Youval Dar
'''

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

    result[file_name] = file_path
  return result

def run():
  '''
  good_MTRIX_pdb_files, good_BIOMT_pdb_files and structure_factors_files
  are dictionaries. the keys are pdb record name and the values are the
  appropriate file full path
  '''
  ## Do the following only if there were changes to the files lists
  #good_MTRIX_pdb_files = make_dict('mtrix_ok_run.txt')
  #good_BIOMT_pdb_files = make_dict('biomt_ok_run.txt')
  #structure_factors_files = make_dict('INDEX','/net/cci/pdb_mirror/structure_factors')
  #pickle.dump(good_MTRIX_pdb_files,open('dict_good_MTRIX_pdb_files','w'))
  #pickle.dump(good_BIOMT_pdb_files,open('dict_good_BIOMT_pdb_files','w'))
  #pickle.dump(structure_factors_files,open('dict_structure_factors_files','w'))

  # If you already have the dictionaries use:
  good_MTRIX_pdb_files = pickle.load(open('dict_good_MTRIX_pdb_files','r'))
  good_BIOMT_pdb_files = pickle.load(open('dict_good_BIOMT_pdb_files','r'))
  structure_factors_files = pickle.load(open('dict_structure_factors_files','r'))
  MTRIX_with_Straucture_Factor = pickle.load(open('MTRIX_with_Straucture_Factor_file_list','r'))

  print 'Dictionaries are loaded...'
  ## When changing the file lists
  #MTRIX_with_Straucture_Factor = []
  i = 0
  for x in good_MTRIX_pdb_files:
    if structure_factors_files.has_key(x):
      i += 1
      MTRIX_with_Straucture_Factor.extend([x,good_MTRIX_pdb_files[x],structure_factors_files[x]])
      #print x
      #print good_MTRIX_pdb_files[x]
      #print structure_factors_files[x]
      #break
  l = len(good_MTRIX_pdb_files)
  print 'The number of both structure factors and good MTRIX with the same name: {} from a total of {}'.format(i,l)
  ## When changing the file lists
  #pickle.dump(MTRIX_with_Straucture_Factor,open('MTRIX_with_Straucture_Factor_file_list','w'))

  i = 0
  for x in good_BIOMT_pdb_files:
    if structure_factors_files.has_key(x):
      i += 1
      #print x
      #print good_BIOMT_pdb_files[x]
      #print structure_factors_files[x]
      #break
  l = len(good_BIOMT_pdb_files)
  print 'The number of both structure factors and good BIOMT with the same name: {} from a total of {}'.format(i,l)

  #i = 0
  #for x in good_MTRIX_pdb_files:
    #print x
    #i += 1
    #if i>5: break
  #print '='*30
  #i = 0
  #for x in structure_factors_files:
    #print x
    #i += 1
    #if i>5: break

  #f1 = open('dict_good_MTRIX_pdb_files.txt','w')
  #f1.writelines([x+'\n' for x in dict_good_MTRIX_pdb_files])
  #f1.close()

  print 'Done...'


if __name__=='__main__':
  #os.chdir('/net/cci-filer2/raid1/home/youval/Work/work')
  os.chdir('C:\Phenix\Dev\Work\work')
  run()
