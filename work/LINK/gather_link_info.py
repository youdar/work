from __future__ import division
from iotbx import pdb
import os

'''
This function collects LINK data from all
PDB files

@author: Youval Dar
'''

def full_pdb_file_paths():
  '''
  Get the list of all PDB files
  '''
  pdb_dir = os.environ["PDB_MIRROR_PDB"]
  pdb_files = open(os.path.join(pdb_dir, "INDEX"), "r").readlines()
  result = []
  for pdb_file in pdb_files:
    pdb_file = os.path.join(pdb_dir, pdb_file.strip())
    # result is a list of records like: '/net/chevy/raid1/pdb_mirror/pdb/00/pdb100d.ent.gz'
    result.append(pdb_file)
  return result

def raw_LINK_data(PDB_files_paths,write_to_file=False):
  '''(list) -> list

  Collect raw info from all PDB files containing LINK info
  add file name to raw data

  Arguments:
  PDB_files_paths: list of all PDB files
  write_to_file: if True, list of file name will be also writen in a file
                 'pdb_file_raw_data.txt'
  '''
  if write_to_file:
    f2 = open('pdb_file_raw_data.txt','w')

  files_with_link_records = []
  for file_name in PDB_files_paths:
    pdb_inp = pdb.input(file_name=file_name)
    link_data = pdb_inp.extract_LINK_records()
    if link_data != []:
      files_with_link_records.append(file_name)
      file_name = file_name.split('/')[-1]
      if write_to_file:
        # add the file name to the end of each data line
        f2.writelines([x+' '+file_name+'\n' for x in link_data])
        f2.writelines('='*30+'\n')

  if write_to_file:
    f2.close()
  return files_with_link_records

def files_with_LINK(files_with_link_records,write_to_file=False):
  '''(list) ->
  '''
  if write_to_file:
    f1 = open('files_with_link_records.txt','w')
    print 'Number of files with LINK records in PDB is: ', len(files_with_link_records)
    f1.writelines([x+'\n' for x in files_with_link_records])
    f1.close()

def run(PDB_files_paths):
  # collect file names of files with LINK info
  files_with_link_records = raw_LINK_data(PDB_files_paths,True)
  # Collect LINK raw data
  files_with_LINK(files_with_link_records)


if __name__=='__main__':
  os.chdir('/net/cci-filer2/raid1/home/youval/Work/work')
  PDB_files_paths = full_pdb_file_paths()
  #pdb_inp = pdb.input(file_name='1W3M.pdb')
  #link_data = pdb_inp.extract_LINK_records()
  run(PDB_files_paths)