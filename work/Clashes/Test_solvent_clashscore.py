from __future__ import division
import mmtbx.monomer_library.pdb_interpretation
from mmtbx import monomer_library
import mmtbx.monomer_library.server
import cctbx.geometry_restraints.manager
from libtbx.utils import Sorry
from libtbx.utils import null_out
#from libtbx.utils import Usage
from libtbx import easy_run
from iotbx import pdb
#import iotbx.utils
import os,sys
import cProfile
import time

''' (str) -> str
Check and collect information on non-bonded clashscore

Argument:
file_name: a pdb file name

Output:
'pdb_file_name::clashscore_all_clashes::clashscore_simple::clashscore_only_sym_op::clashscore_solvent_solvent'

>>> python Test_solvent_clashscore.py 1a18.pdb
1a18::13.7::10.1::10.0310::10.1070
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
  pdb_processed_file = monomer_library.pdb_interpretation.run(args=[file_name],
    substitute_non_crystallographic_unit_cell_if_necessary=True,
    assume_hydrogens_all_missing=False,
    hard_minimum_nonbonded_distance=0.0,
    nonbonded_distance_threshold=None,
    log=null_out())
  
  pdb_inp = pdb.hierarchy.input(file_name=file_name)
  
  grm = pdb_processed_file.geometry_restraints_manager(assume_hydrogens_all_missing=False)

  xrs = pdb_processed_file.xray_structure()
  sites_cart = xrs.sites_cart()
  site_labels = xrs.scatterers().extract_labels()
  hd_sel = xrs.hd_selection()

  obj = grm.get_nonbonded_clashscore(
    sites_cart=sites_cart,
    site_labels=site_labels,
    hd_sel=hd_sel)

  clashscore_all = obj.nb_clashscore_all_clashes
  clashscore_simple = obj.nb_clashscore_simple
  clashscore_due_to_sym_op = obj.nb_clashscore_due_to_sym_op
  clashscore_solvent_solvent = obj.nb_clashscore_solvent_solvent

  return obj

def get_new_file_name(file_name):
  '''(str) -> str
  When using ready_set or reduce, a new file will be created, one that includes
  hydrogens. This function finds the new name and returns it

  Argument:
  file_name: a file name of a pdb file, including the path, that was processed
  using ready_set or reduce

  Return:
  file_name: the file name of the modified pdb file, the one that includes hydrogens
  '''
  tmp = file_name.split('.')
  if not (len(tmp)==3 and tmp[1] == 'updated'):
    file_name = get_pdb_file(file_name,fetch_file=True,print_out=False)
    # when using phenix.reduce
    cmd = 'phenix.reduce {0} > {1}.updated.{2}'.format(file_name,tmp[0],tmp[-1])
    
    # when using ready_set
    #cmd  = 'phenix.ready_set {0}'.format(file_name)
    #
    probe_out = easy_run.go(cmd)
    # getting to modified file name
    # when using phenix.reduce
    file_name = '{0}.updated.{1}'.format(tmp[0],tmp[-1])
    
    # when using ready_set
    #line  = [x for x in probe_out.stdout_lines if 'Writing model to' in x]
    #if line != []:
      #file_name = line[0].split('Writing model to ')[-1]
    #else:
      #collect_sorry = [x for x in probe_out.stdout_lines if 'Sorry:' in x]
      #for x in collect_sorry:
        #print x

  return file_name

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

def get_pdb_file(file_name,fetch_file=False,print_out=True):
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
    if print_out:
      print 'No such file in working directory. Trying to fetch {} file from RCSB web site'.format(file_name)
    file_name = get_file_name(file_name)
    osType = sys.platform
    if osType.startswith('win') or fetch_file==True:
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
  elif print_out:
    print 'Using the file {} found in the working directory'.format(file_name)
  return file_name

def run(file_name):
  set_working_path('Clashes\wtest')
  file_name = sys.argv[1]
  # get the file name of the file that includes hydrogens
  file_name = get_new_file_name(file_name)
  #
  obj = get_clashscore_internal(file_name)
  # Cleanup
  #if os.path.isfile(file_name):
    #os.remove(file_name)
  #file_name_2 = file_name[0:5] + file_name[-3:]
  #if os.path.isfile(file_name_2):
    #os.remove(file_name_2)
  #print '\nClashscore all: {0:.2f}\nsimple: {1:.2f}\nonly_sym_op   : {2:.3f}\n'.format(*clashscore)
  output_file_name = get_file_name(file_name)
  outstr = '{0}::{1:.1f}::{2:.1f}::{3:.1f}::{4:.1f}'.format(
    output_file_name,  					# pdb 4 letter file name
    obj.nb_clashscore_all_clashes,		# clashscore_all_clashes
    obj.nb_clashscore_simple,  			# clashscore_simple
    obj.nb_clashscore_due_to_sym_op,	# clashscore_only_sym_op
    obj.nb_clashscore_solvent_solvent) # clashscore_solvent_solvent
  print outstr
  
if (__name__ == "__main__"):
  file_name = sys.argv[1]
  run(file_name)




