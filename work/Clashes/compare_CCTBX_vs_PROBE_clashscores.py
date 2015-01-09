from __future__ import division
from libtbx.utils import null_out
from libtbx.utils import Sorry
from libtbx import easy_run
from iotbx.pdb import fetch
import getpass
import os,sys


def get_pdb_id(file_name):
  '''
  clean a pdb file name, remove path and file extensions

  Args:
    file_name (str): file name that might include path and extension

  Returns:
    file_name (str): 4 letter pdb id
  '''
  file_name = os.path.split(file_name)[-1]
  file_name = file_name.lower()
  file_name = file_name.split('.')[0]
  if len(file_name)>4:
    if 'pdb' in file_name:
      i = file_name.find('pdb')
      file_name = file_name[i+3:i+7]
  return file_name

def get_file_from_rcsb(pdb_id,format='pdb'):
  """ (file_name) -> file_path
  fetch pdb or structure factor file for pdb_id from the RCSB website

  Args:
    file_name: a pdb file name
    format (str): for PDB 'pdb' for CIF 'cif'

  Returns:
    a file path for the pdb file_name
  """
  try:
    file_name = fetch.get_pdb(
      pdb_id,
      data_type='pdb',
      format=format,
      mirror='rcsb',
      log=null_out())
  except Sorry:
    file_name = ''
  return file_name

def run(file_name,verbose=False):
  """
  Update file_name to include hydrogens, then get clashscores for CCTBX and
  PROBE

  Args:
    file_name (str): pdb file name including path that should be
      processed using ready_set or reduce
    verbose (bool): when True, return results in a table

  Returns:
    out_str (str):
      if verbose is False
        x1::x2::x3::x4::x5::x6
          Where:
            x1: pdb_file_name
            x2: simple_cctbx_clashscore
            x3: symmetry_cctbx_clashscore
            x4: solvent_cctbx_clashscore
            x5: total_cctbx_clashscore
            x6: probe_clashscore

      if verbose is True
             PDB id  |        CCTBX                     | PROBE
                     | Simple   Sym      Solvent  All   |
            ----------------------------------------------------
              x1     |  x2       x3       x4      x5    | x6
  """
  _,fn = os.path.split(file_name)
  pdb_id = get_pdb_id(fn)
  # Get files in pdb format even when at LBL, to avoid issues with phenix.reduce
  if ('_mirror' in file_name) or (file_name == pdb_id):
    file_name = pdb_id + '.pdb'
  file_in_local_folder = os.path.isfile(file_name)
  file_to_clean = []
  if not file_in_local_folder:
    # leave file in folder if it was already there
    file_name = get_file_from_rcsb(pdb_id,format='pdb')
    file_to_clean.append(file_name)
  file_with_hydrogens = file_name.replace('.pdb','_with_h.pdb')
  file_to_clean.append(file_with_hydrogens)
  # get cif file
  cif_file = get_file_from_rcsb(pdb_id,format='cif')
  file_to_clean.append(cif_file)
  #
  if os.path.isfile(file_name):
    cmd = 'phenix.clashscore {} verbose=false method={}'
    files = file_name + ' ' + cif_file
    r = easy_run.go(cmd.format(file_name,'molprobity'),join_stdout_stderr=False)
    r2 = easy_run.go(cmd.format(files,'cctbx'),join_stdout_stderr=False)

    if bool(r.stderr_lines) or bool(r2.stderr_lines):
      # Error during processing
      simple = sym = solv = all = clashscore_probe = -2
    else:
      clashscore_probe = float(r.stdout_lines[0])
      # clean cctbx results
      r2 = r2.stdout_lines[0]
      r2 = r2.replace('[','')
      r2 = r2.replace(']','')
      simple,sym,solv,all = [float(x) for x in r2.split(',')]
  else:
    # file_name is not a valid pdb file
    simple = sym = solv = all = clashscore_probe = -1
  # Cleanup files if they where not already in local folder
  if file_to_clean:
    for fn in file_to_clean:
      if os.path.isfile(fn): os.remove(fn)
  #
  pdb_id = get_pdb_id(file_name)
  data = [
    pdb_id,                     # pdb 4 letter file name
    round(simple,1),  	        # CCTBX clashscore_simple
    round(sym,1),	              # CCTBX clashscore_only_sym_op
    round(solv,1),              # CCTBX clashscore_solvent_solvent
    round(all,1),               # CCTBX clashscore_all_clashes
    round(clashscore_probe,1)]  # PROBE clashscore
  if verbose:
    table_str = '{:^8} | {:<7}  {:<7}  {:<7}  {:<5} | {:<5}'
    title1 = ['PDB id','','CCTBX','','','PROBE']
    title2 = ['','Simple','Sym','Solvent','All','']
    outstr = '\n' + table_str.format(*title1) + '\n'
    outstr += table_str.format(*title2) + '\n'
    outstr += '-'*52 + '\n'
    outstr += table_str.format(*data) + '\n'
  else:
    outstr = '{0}::{1:.1f}::{2:.1f}::{3:.1f}::{4:.1f}::{5:.1f}'.format(*data)
  return outstr

def set_working_folder():
  """ set working folder path for Youval """
  username = getpass.getuser()
  osType = sys.platform
  if username.lower() == 'youval':
    if osType.startswith('win'):
      dr = r'C:\Phenix\Dev\Work\work\Clashes\wtest'
    else:
      dr = '/net/cci-filer2/raid1/home/youval/Work/work/Clashes/wtest'
    os.chdir(dr)

if (__name__ == "__main__"):
  help_str = """\
Get clashscores for CCTBX and PROBE after updating the pdb file to include
hydrogens

Examples:
>>>python compare_clashscores_CCTBX_vs_PROBE.py 1a18.pdb verbose
 PDB id  |          CCTBX                   | PROBE
         | Simple   Sym      Solvent  All   |
----------------------------------------------------
  1a18   | 13.7     2.4      0.0      16.1  | 13.0

>>>python compare_clashscores_CCTBX_vs_PROBE.py 1a18.pdb
1a18::16.1::13.7::2.4::0.0::13.0
  """
  args = sys.argv[1:]
  if len(args) == 0:
    print help_str
  else:
    current_folder = os.getcwd()
    file_name = args[0]
    verbose = False
    if len(args) > 1:
      verbose = [('verbose' in x) or ('-v' == x) for x in args[1:]]
      verbose = (verbose.count(True) > 0)
    set_working_folder()
    print run(file_name,verbose=verbose)
    os.chdir(current_folder)






