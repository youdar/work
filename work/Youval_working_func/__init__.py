from __future__ import division
import sys
import os

class null_out(object):
  """Pseudo-filehandle for suppressing printed output."""
  def isatty(self): return False
  def close(self): pass
  def flush(self): pass
  def write(self, str): pass
  def writelines(self, sequence): pass


def set_working_path(work_path):
  """
  Change testing directory to work_path

  Note:  the base working folders are specific to my setup

  >>>set_working_path('Clashes\wtest')
  """
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
    directory_path = '~/youval/Work/work/{0}'.format(work_path)
  os.chdir(directory_path)


def get_4_letters_pdb_id(file_name):
  """(str)  -> str
  clean a pdb file name, remove path and file extensions

  :param file_name (str): pdb file name that may look like pdb1a37.pdb
  :return pdb_id (str): the 4 letter pdb id

  >>>get_4_letters_pdb_id('pdb1a37.pdb')
  1a37
  >>>get_4_letters_pdb_id('1a37')
  1a37
  """
  basename = os.path.basename(file_name)
  file_name, file_type = os.path.splitext(basename)
  if len(file_name)>4:
    if 'pdb' in file_name:
      i = file_name.find('pdb')
      pdb_id = file_name[i+3:i+7]
  elif len(file_name)==4:
    pdb_id = file_name
  else:
    pdb_id = None
  return pdb_id


def get_pdb_file(file_name, print_out=True):
  """ (file_name) -> file_path
  This function will check if a pdb file_name exist.
  If it is not, it will either fetch it from the RCSB website
  or find it on LBLs pdb mirror folder

  :param file_name (str): a pdb file name
  :return file_name (str): the location, path, of the file
  """
  if not os.path.isfile(file_name):
    # get a clean pdb file name
    if print_out:
      s = 'No such file in working directory. ' \
          'Trying to fetch {} file from RCSB web site'
      print s.format(file_name)
    file_name = get_4_letters_pdb_id(file_name)
    osType = sys.platform
    if osType.startswith('win'):
      from iotbx.pdb import fetch
      # fetch pdb file from intenet
      file_name = fetch.get_pdb (file_name,'pdb',mirror='rcsb',log=null_out())
    else:
      # find the file in LBL pdb mirror folder
      pdb_dir = os.environ["PDB_MIRROR_PDB"]
      pdb_files = open(os.path.join(pdb_dir, "INDEX"), "r").read().splitlines()
      for i,p in enumerate(pdb_files):
        if file_name in p:
          break
      file_name = os.path.join(pdb_dir,pdb_files[i])
      # # Second method
      # f = 'pdb{}.ent.gz'.format(file_name)
      # file_name = []
      # for root, _, files in os.walk(pdb_dir):
      #   if f in files:
      #     file_name = os.path.join(root, f)
      #     break
  elif print_out:
    print 'Using the file {} found in the working directory'.format(file_name)
  return file_name
