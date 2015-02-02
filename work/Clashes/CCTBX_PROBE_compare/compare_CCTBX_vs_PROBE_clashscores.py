from __future__ import division
import mmtbx.monomer_library.pdb_interpretation as pdb_inter
import cctbx.geometry_restraints.nonbonded_overlaps as cs
from libtbx.utils import null_out
from libtbx import easy_run
import iotbx.pdb
import getpass
import sys
import os


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

def run(file_name,verbose=False,macro_mol_mode=False):
  """
  Update file_name to include hydrogens, then get clashscores for CCTBX and
  PROBE

  Args:
    file_name (str): pdb file name including path that should be
      processed using ready_set or reduce
    verbose (bool): when True, return results in a table
    macro_mol_mode (bool): when True, clean PDB file, leaving only protein,
      RNA or DNA.

  Returns:
    out_str (str):
      if verbose is False
        x1::x2::x3::x4::x5::x6
          Where:
            x1: pdb_file_name
            x2: macro_molecule_cctbx_clashscore
            x3: symmetry_cctbx_clashscore
            x4: total_cctbx_clashscore
            x5: probe_clashscore

      if verbose is True
              PDB id  |             CCTBX          | PROBE
                      | Macro Mol.  Sym      All   |
              ----------------------------------------------
               x1     | x2          x3       x4    | x5

    Error Types:
      -2: Could process CCTBX but not PROBE
      -3: Could process PROBE but not CCTBX (model contains unknown_type_pairs)
      -4: Could process PROBE but not CCTBX (multiple models in pdb files)
      -5: Could process PROBE but not CCTBX (Other processing issue)
      -6: Bad CRYST1 records, bad crystal symmetry
      -7: Both CCTBX and PROBE fail processing
  """
  _,fn = os.path.split(file_name)
  pdb_id = get_pdb_id(fn)
  cif_file = ''
  # Get files in pdb format even when at LBL, to avoid issues with phenix.reduce
  if ('_mirror' in file_name) or (file_name == pdb_id):
    file_name = pdb_id + '.pdb'
  file_to_clean = []
  if not os.path.isfile(file_name):
    # leave file in folder if it was already there
    # cmd = 'phenix.fetch_pdb {} --all'.format(pdb_id)
    cmd = 'phenix.fetch_pdb {}'.format(pdb_id)
    r = easy_run.go(cmd,join_stdout_stderr=False)
    for fn in r.stdout_lines:
      fn = os.path.split(fn)[-1]
      # if '.cif' in fn: cif_file = fn
      if '.pdb' in fn: file_name = fn
      file_to_clean.append(fn)
  if macro_mol_mode:
    file_name = save_macro_mol_mode(file_name)
    file_to_clean.append(file_name)
  file_with_hydrogens = file_name.replace('.pdb','_with_h.pdb')
  file_to_clean.append(file_with_hydrogens)
  #
  if os.path.isfile(file_name):
    option = 'substitute_non_crystallographic_unit_cell_if_necessary=false'
    cmd = 'phenix.clashscore {} verbose=false method={} ' + option
    files = file_name + ' ' + cif_file
    r = easy_run.go(cmd.format(file_name,'molprobity'),join_stdout_stderr=False)
    r2 = easy_run.go(cmd.format(files,'cctbx'),join_stdout_stderr=False)
    cctbx_error = bool(r2.stderr_lines)
    probe_error = bool(r.stderr_lines)
    if not cctbx_error and probe_error:
      # Error processing molprobity
      macro_molecule = sym = all = clashscore_probe = -2
    elif not probe_error and cctbx_error:
      # Error processing cctbx
      unknown_type_pairs = False
      multiple_models = False
      bad_cryst1 = False
      for l in r2.stderr_lines:
        unknown_type_pairs |= 'unknown type pairs' in l
        multiple_models |= 'provide only a single model' in l
        bad_cryst1 |= 'None valid CRSYT1 records' in l
        if unknown_type_pairs:
          macro_molecule = sym = all = clashscore_probe = -3
        elif multiple_models:
          macro_molecule = sym = all = clashscore_probe = -4
        elif bad_cryst1:
          macro_molecule = sym = all = clashscore_probe = -5
        else:
          macro_molecule = sym = all = clashscore_probe = -6
    elif probe_error and cctbx_error:
      # error on both PROBE and CCTBX
      macro_molecule = sym = all = clashscore_probe = -7
    else:
      clashscore_probe = float(r.stdout_lines[0])
      # clean cctbx results
      r2 = r2.stdout_lines[0]
      r2 = r2.replace('[','')
      r2 = r2.replace(']','')
      macro_molecule,sym,all = [float(x) for x in r2.split(',')]
  else:
    # file_name is not a valid pdb file
    macro_molecule = sym = all = clashscore_probe = -1
  # Cleanup files if they where not already in local folder
  # if file_to_clean:
  #   for fn in file_to_clean:
  #     if os.path.isfile(fn): os.remove(fn)
  #
  pdb_id = get_pdb_id(file_name)
  data = [
    pdb_id,                     # pdb 4 letter file name
    round(macro_molecule,1),  	# CCTBX clashscore_macro_molecule
    round(sym,1),	              # CCTBX clashscore_only_sym_op
    round(all,1),               # CCTBX clashscore_all
    round(clashscore_probe,1)]  # PROBE clashscore
  if verbose:
    error_dict = {
      -2: 'Could process CCTBX but not PROBE',
      -3:
        'Could process PROBE but not CCTBX (model contains unknown_type_pairs)',
      -4: 'Could process PROBE but not CCTBX (multiple models in pdb files)',
      -5: 'Could process PROBE but not CCTBX (Other processing issue)',
      -6: 'Both CCTBX and PROBE fail processing'}
    table_str = '{:^8} | {:<10}  {:<7}  {:<5} | {:<5}'
    title1 = ['PDB id','','CCTBX','','PROBE']
    title2 = ['','Macro Mol.','Sym','All','']
    outstr = '\n' + table_str.format(*title1) + '\n'
    outstr += table_str.format(*title2) + '\n'
    outstr += '-'*46 + '\n'
    if macro_molecule < 0:
      outstr += error_dict[macro_molecule]
    else:
      outstr += table_str.format(*data) + '\n'
  else:
    outstr = '{0}::{1:.1f}::{2:.1f}::{3:.1f}::{4:.1f}'.format(*data)
  return outstr

def save_macro_mol_mode(file_name):
  """
  Create a new pdb file that contains only the macro molecule

  Args:
    file_name (str): PDB file name
  """
  assert os.path.isfile(file_name)
  fn = file_name.replace('.pdb','_macro_mol.pdb')
  macro_molecule_selection_str = 'protein or dna or rna'

  pdb_processed_file = pdb_inter.run(
    args=[file_name],
    assume_hydrogens_all_missing=False,
    hard_minimum_nonbonded_distance=0.0,
    nonbonded_distance_threshold=None,
    substitute_non_crystallographic_unit_cell_if_necessary=True,
    log=null_out())
  macro_mol_sel = cs.get_macro_mol_sel(pdb_processed_file)
  pdb_inp = iotbx.pdb.input(file_name=file_name)
  ph = pdb_inp.construct_hierarchy()
  ph = ph.select(macro_mol_sel)
  cryst_sym = pdb_inp.crystal_symmetry()
  pdb_str = ph.as_pdb_string(cryst_sym)
  open(fn,'w').write(pdb_str)
  return fn

def set_working_folder():
  """ set working folder path for Youval """
  username = getpass.getuser()
  osType = sys.platform
  if username.lower() == 'youval':
    if osType.startswith('win'):
      dr = r'C:\Phenix\Dev\Work\work\Clashes\wtest'
    else:
      dr = '/net/cci/youval/work/work/Clashes/wtest'
    os.chdir(dr)

if (__name__ == "__main__"):
  help_str = """\
Get clashscores for CCTBX and PROBE after updating the pdb file to include
hydrogens

Examples:
>>>python compare_clashscores_CCTBX_vs_PROBE.py 1a18.pdb verbose
 PDB id  |             CCTBX          | PROBE
         | Macro Mol.  Sym      All   |
----------------------------------------------
  1a18   | 13.2        0.0      13.7  | 13.0

>>>python compare_clashscores_CCTBX_vs_PROBE.py 1a18.pdb
1a18::16.1::13.7::2.4::13.0
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
      macro_mol_mode = [('clean' in x) or ('-c' == x) for x in args[1:]]
      macro_mol_mode = (macro_mol_mode.count(True) > 0)
    set_working_folder()
    print run(file_name,verbose=verbose,macro_mol_mode=macro_mol_mode)
    os.chdir(current_folder)






