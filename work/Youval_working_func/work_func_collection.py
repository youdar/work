from __future__ import division
import os, sys

'''
A collection of helping sunction I use regularly 

@Author: Youval Dar
'''

def set_working_path(work_path):
  '''(str)
  Change working directory to work_path
  
  Note:  the base working folders are specific to my setup
  
  >>>set_working_path('Clashes\junk')
  '''
  # locate the directory containing the log files
  if '\\' in work_path:
    work_path = work_path.split('\\')
  elif '/' in work_path:
    work_path = work_path.split('/')
  osType = sys.platform
  if osType.startswith('win'):
    work_path = '\\'.join(work_path)
    directory_path = 'c:\\Phenix\\Dev\\Work\\work\\{0}'.format(work_path)
  else:
    work_path = '/'.join(work_path)
    directory_path = '/net/cci-filer2/raid1/home/youval/Work/work/{0}'.format(work_path)
  os.chdir(directory_path)
  
def get_file_name(file_name):
  '''(str)  -> str
  clean a pdb file name, remove path and file extensions
  
  Return the 4 letter pdb id
  '''
  osType = sys.platform
  if osType.startswith('win'):
      file_name = file_name.split('\\')[-1]
  else:
      file_name = file_name.split('/')[-1]
  file_name = file_name.split('.')[0]
  if len(file_name)>4:
    if 'pdb' in file_name:
      i = file_name.find('pdb')
      file_name = file_name[i+3:i+7]
  return file_name


def get_pdb_file(file_name):
  ''' (file_name) -> file_path
  This function will check if a pdb file_name exist.
  If it is not, it will either fetch it from the RCSB website 
  or find it on LBLs pdb mirror folder

  Argument:
  file_name: a pdb file name
  
  Return:
  a file path for the pdb file_name
  '''
  from iotbx.pdb import fetch
  class null_out(object):
    """Pseudo-filehandle for suppressing printed output."""
    def isatty(self): return False
    def close(self): pass
    def flush(self): pass
    def write(self, str): pass
    def writelines(self, sequence): pass
  log  = null_out()
  
  if not os.path.isfile(file_name):
    # get a clean pdb file name
    print 'No such file in working directory. Try to fetch {} file from RCSB web site'.format(file_name)
    file_name = get_file_name(file_name)  
    osType = sys.platform
    if osType.startswith('win'):
      # fetch pdb file
      file_name = fetch.get_pdb (file_name,'pdb',mirror='rcsb',log=log)
    else:
      # find the file in LBL pdb mirror folder
      pdb_dir = os.environ["PDB_MIRROR_PDB"]
      pdb_files = open(os.path.join(pdb_dir, "INDEX"), "r").read().splitlines()
      for i,p in enumerate(pdb_files):
        if file_name in p:
          break
      file_name = os.path.join(pdb_dir,pdb_files[i])
  else:
    print 'Using the file {} found in the working directory'.format(file_name)
  return file_name