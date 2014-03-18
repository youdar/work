from __future__ import division
import cPickle as pickle
from libtbx import smart_open
import os,sys
import datetime



def run():
  data_file_path = '/net/cci/youval/Work/work/Clashes/Data'
  file_to_year_dict = {}

  # find the file in LBL pdb mirror folder
  #pdb_dir = os.environ["PDB_MIRROR_PDB"]
  #pdb_files = open(os.path.join(pdb_dir, "INDEX"), "r").read().splitlines()
  #for p in pdb_files:
    #file_name = p[-11:-7]
    #file_name_with_path = os.path.join(pdb_dir,p)
    #file_lines = smart_open.for_reading(
            #file_name = file_name_with_path).read().splitlines()
    #year = get_year(file_lines)
    #file_to_year_dict[file_name] = year

  #print len(file_to_year_dict)

  # collect year for files that are not in INDEX
  os.chdir('/net/cci/youval/Work/work/Clashes/junk')
  from iotbx.pdb import fetch
  class null_out(object):
    """Pseudo-filehandle for suppressing printed output."""
    def isatty(self): return False
    def close(self): pass
    def flush(self): pass
    def write(self, str): pass
    def writelines(self, sequence): pass
  missing_files = pickle.load(open(os.path.join(data_file_path,'missing_files'),'r'))
  for f in missing_files:
    file_name = fetch.get_pdb (f,'pdb',mirror='rcsb',log=null_out())
    file_lines = smart_open.for_reading(
                file_name = file_name).read().splitlines()
    year = get_year(file_lines)
    file_to_year_dict[f] = year
    print len(file_to_year_dict)
    os.remove(file_name)

  #pickle.dump(file_to_year_dict, open(os.path.join(data_file_path, "file_to_year_dict"),'w'))
  print len(file_to_year_dict)
  print 'Done...'

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