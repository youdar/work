from __future__ import division
import cPickle as pickle
from libtbx import smart_open
import os,sys
import datetime



def run():
  data_dir = '/net/cci/youval/Work/work/MTRIX/Data'
  #data_dir = r'c:\Phenix\Dev\Work\work\MTRIX\Data'
  os.chdir(data_dir)
  file_to_year_dict = {}

  files_with_good_MTRIX = set(pickle.load(open(os.path.join(data_dir,'files_with_good_MTRIX'),'r')))
  good_MTRIX_pdb_files = pickle.load(open(os.path.join(data_dir,'dict_good_MTRIX_pdb_files'),'r'))

  # find the file in LBL pdb mirror folder
  for fn in files_with_good_MTRIX:
    file_name_with_path = good_MTRIX_pdb_files[fn]
    file_lines = smart_open.for_reading(
            file_name = file_name_with_path).read().splitlines()
    year = get_year(file_lines)
    file_to_year_dict[fn] = year

  print len(file_to_year_dict)

  ##pickle.dump(file_to_year_dict, open(os.path.join(data_dir, "file_to_year_dict"),'w'))
  #print len(file_to_year_dict)
  #print 'Done...'

def get_year(file_lines):
  for line in file_lines:
    if line.startswith('HEADER'):
      year = line[57:59]
      if year > '30':
        year = '19' + year
      else:
        year = '20' + year
      return year
      break
  return ''


if (__name__ == "__main__"):
  run()