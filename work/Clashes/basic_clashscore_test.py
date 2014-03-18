from __future__ import division
import mmtbx.monomer_library.pdb_interpretation
from mmtbx import monomer_library
import mmtbx.validation.clashscore
import mmtbx.monomer_library.server
import cctbx.geometry_restraints.manager
from libtbx.utils import Sorry
from libtbx.utils import null_out
from libtbx.utils import Usage
from iotbx import pdb
import iotbx.utils
import iotbx.phil
import os,sys
import cProfile
import time

''' (str) -> str
Compare clash scores of phenix.clashscore and the internal clashscore function

Argument:
file_name: a pdb file name

Output:
'pdb_file_name::total_nb_clashscore:::without_sym_nb_clashscore::clashscore_probe'

>>> python Test_internal_clashscore.py 1a18.pdb [10,11,'NMR'] 1a18
1a18::13.7::10.1::10
'''

def get_clashscore_internal(file_name):
  '''(str) -> float

  Calculate clashscore using pdb_interpertation code.

  Argument:
  file_name: a string, a file name, including the path, of a pdb file
  to which hydorgen were added

  Output:
  clashscore: a float number representning the clashscore of the pdb file

  '''
  pdb = monomer_library.pdb_interpretation.run(
    args=[file_name],
    substitute_non_crystallographic_unit_cell_if_necessary=True,
    assume_hydrogens_all_missing=False,
    hard_minimum_nonbonded_distance=0.0,
    log=null_out())

  grm = pdb.geometry_restraints_manager(
    assume_hydrogens_all_missing=False,
    hard_minimum_nonbonded_distance=0.0)

  xrs = pdb.xray_structure()
  sites_cart = xrs.sites_cart()
  site_labels = xrs.scatterers().extract_labels()
  hd_sel = xrs.hd_selection()

  obj = grm.get_nonbonded_clashscore(
    sites_cart=sites_cart,
    site_labels=site_labels,
    hd_sel=hd_sel)

  clashscore_all = obj.nb_clashscore_all_clashes
  clashscore_without_sym_op = obj.nb_clashscore_without_sym_op
  clashscore_due_to_sym_op = obj.nb_clashscore_due_to_sym_op

  return clashscore_all,clashscore_without_sym_op,clashscore_due_to_sym_op


def set_working_path(work_path):
  '''(str)
  Change working directory to work_path

  Note:  the base working folders are specific to my setup

  >>>set_working_path('Clashes\wtest')
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

#def get_pdb_file(file_name, print_out=True):
  #''' (file_name) -> file_path
  #This function will check if a pdb file_name exist.
  #If it is not, it will either fetch it from the RCSB website
  #or find it on LBLs pdb mirror folder

  #Argument:
  #file_name: a pdb file name

  #Return:
  #a file path for the pdb file_name
  #'''
  #from iotbx.pdb import fetch
  #class null_out(object):
    #"""Pseudo-filehandle for suppressing printed output."""
    #def isatty(self): return False
    #def close(self): pass
    #def flush(self): pass
    #def write(self, str): pass
    #def writelines(self, sequence): pass
  #log  = null_out()

  #if not os.path.isfile(file_name):
    ## get a clean pdb file name
    #if print_out:
      #print 'No such file in working directory. Trying to fetch {} file from RCSB web site'.format(file_name)
    #file_name = get_file_name(file_name)
    #osType = sys.platform
    #if osType.startswith('win'):
      ## fetch pdb file
      #file_name = fetch.get_pdb (file_name,'pdb',mirror='rcsb',log=log)
    #else:
      ## find the file in LBL pdb mirror folder
      #pdb_dir = os.environ["PDB_MIRROR_PDB"]
      #pdb_files = open(os.path.join(pdb_dir, "INDEX"), "r").read().splitlines()
      #for i,p in enumerate(pdb_files):
        #if file_name in p:
          #break
      #file_name = os.path.join(pdb_dir,pdb_files[i])
  #elif print_out:
    #print 'Using the file {} found in the working directory'.format(file_name)
  #return file_name


if (__name__ == "__main__"):
  set_working_path('Clashes\junk')
  file_name = sys.argv[1]
  clashscore_all,clashscore_without_sym_op,clashscore_only_sym_op = get_clashscore_internal(file_name)
  # probe clashscore
  score_with_hydrogen,score_without_hydrogen,experment_type = sys.argv[2]
  # score_without_hydrogen means that keep_hydrogens=False

  outstr = '{0}::{1:.1f}::{2:.1f}::{3:.1f}'.format(
    sys.argv[3],
    clashscore_all,
    clashscore_without_sym_op,
    score_with_hydrogen)
  print outstr




