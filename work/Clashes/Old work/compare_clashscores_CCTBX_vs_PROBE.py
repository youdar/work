from __future__ import division
import mmtbx.monomer_library.pdb_interpretation
import mmtbx.monomer_library.server
import mmtbx.validation.clashscore
from mmtbx import monomer_library
from libtbx.utils import null_out
from libtbx.utils import Sorry
from libtbx.utils import Usage
from libtbx import easy_run
from iotbx.pdb import fetch
import iotbx.utils
import iotbx.phil
import getpass
import os,sys


def get_cctbx_clashscores(file_name, cif_file=None):
  """
  Calculate clashscore using pdb_interpertation code.

  Args:
  file_name (str): pdb file name including path, to which hydrogens were added
    (adding hydrogens is done using the get_updated_file_name function)

  Returns:
    simple (float): clashscore not due to symmetry or solvent-solvent clashes
    sym (float): clashscore due to symmetry operations
    solv (float): clashscore due to solvent-solvent clashes
    all (float): sum of all three clashscores above
  """
  args = [file_name]
  if cif_file:
    args.append(cif_file)
  pdb = monomer_library.pdb_interpretation.run(
    args=args,
    substitute_non_crystallographic_unit_cell_if_necessary=True,
    assume_hydrogens_all_missing=False,
    hard_minimum_nonbonded_distance=0.0,
    nonbonded_distance_threshold=None,
    log=null_out())
  # build restraints manager and x-ray structure
  grm = pdb.geometry_restraints_manager(assume_hydrogens_all_missing=False)
  xrs = pdb.xray_structure()
  sites_cart = xrs.sites_cart()
  site_labels = xrs.scatterers().extract_labels()
  hd_sel = xrs.hd_selection()
  # get cctbx non-bonded clashscore object
  try:
    obj = grm.get_nonbonded_clashscore(
      sites_cart=sites_cart,
      site_labels=site_labels,
      hd_sel=hd_sel)
    #
    all = obj.nb_clashscore_all_clashes
    simple = obj.nb_clashscore_simple
    sym = obj.nb_clashscore_due_to_sym_op
    solv = obj.nb_clashscore_solvent_solvent
  except Sorry:
    simple = sym = solv = all = -2
  return simple,sym,solv,all

def get_updated_file_name(file_name,out_folder='',use_reduce=True):
  """
  When using ready_set or reduce, a new file will be created, one that includes
  hydrogens. This function finds the new name and returns it

  Args:
    file_name: pdb file name including path that should be
      processed using ready_set or reduce
    use_reduce (bool): When True use phenix.reduce, otherwise phenix.ready_set
    out_folder (str): the folder where the updated file will be writen (when
      different than the source file location)

  Return:
    file_name: the file name of the modified pdb file, the one that includes
      hydrogens
  """
  (path,file_name) = os.path.split(file_name)
  pdb_id = get_pdb_id(file_name)
  if out_folder:
    path = out_folder
  tmp = file_name.split('.')
  # check that the file is not already an updated file
  if not (len(tmp)==3 and tmp[1] == 'updated'):
    file_name = get_file_from_rcsb(pdb_id)
    # stop processing if file_name was not found
    if not os.path.isfile(file_name): return ''
    reduce = 'phenix.reduce {0} > {1}.updated.pdb'
    ready_set = '{0}.updated.{1}'
    if use_reduce:
      #file_name = file_name[-8:]
      cmd = reduce.format(file_name,pdb_id)
    else:
      cmd  = 'phenix.ready_set {0}'.format(file_name)
    #
    probe_out = easy_run.go(cmd)
    # getting to modified file name
    if use_reduce:
      file_name = ready_set.format(tmp[0],tmp[-1])
    else:
      line  = [x for x in probe_out.stdout_lines if 'Writing model to' in x]
      if line != []:
        file_name = line[0].split('Writing model to ')[-1]
      else:
        collect_sorry = [x for x in probe_out.stdout_lines if 'Sorry:' in x]
        for x in collect_sorry:
          print x
  return file_name

def get_probe_clashscore(file_name):
  """
  Calculate PROBE clashscore.

  Args:
   file_name (str): pdb file name including path, to which hydrogens were added
   (adding hydrogens is done using the get_updated_file_name function)

  Returns:
    clashscore (float):PROBE clashscore
  """
  usage_string ='''\
phenix.clashscore file.pdb [params.eff] [options ...]

Options:

  model=input_file      input PDB file
  keep_hydrogens=True   keep input hydrogen files (otherwise regenerate)
  nuclear=False         use nuclear x-H distances and vdW radii
  verbose=True          verbose text output

Example:

  phenix.clashscore model=1ubq.pdb keep_hydrogens=True

'''
  cmdline = iotbx.phil.process_command_line_with_files(
    args=[file_name],
    master_phil=get_master_phil(),
    pdb_file_def="model",
    usage_string=usage_string)
  params = cmdline.work.extract()
  if (params.model is None) :
    raise Usage(usage_string)
  pdb_in = cmdline.get_file(params.model, force_type="pdb")
  hierarchy = pdb_in.file_object.construct_hierarchy()
  # We keep hydrogen
  result_with_hydrogen = mmtbx.validation.clashscore.clashscore(
    pdb_hierarchy=hierarchy,
    keep_hydrogens=True,
    nuclear=params.nuclear,
    out=null_out(),
    verbose=False)
  return result_with_hydrogen.get_clashscore()

def get_master_phil() :
  # phil parameters for the probe clashscore
  phil_str = """\
model = None
  .type = path

verbose = True
  .type = bool

keep_hydrogens = True
  .type = bool
  .help = '''Keep hydrogens in input file'''

nuclear = False
  .type = bool
  .help = '''Use nuclear hydrogen positions'''

time_limit = 120
  .type = int
  .help = '''Time limit (sec) for Reduce optimization'''
  """
  return iotbx.phil.parse(phil_str)

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

def get_file_from_rcsb(pdb_id,data_type='pdb'):
  """ (file_name) -> file_path
  fetch pdb or structure factor file for pdb_id from the RCSB website

  Args:
    file_name: a pdb file name
    data_type (str):
      'pdb' -> pdb
      'xray' -> structure factor

  Returns:
    a file path for the pdb file_name
  """
  try:
    file_name = fetch.get_pdb(pdb_id,data_type,mirror='rcsb',log=null_out())
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
    file_to_clean.append(file_name)
  file_has_hydrogens = ('.updated.' in file_name)
  # get the file name of the file that includes hydrogens
  if not file_has_hydrogens:
    file_name = get_updated_file_name(file_name)
    file_to_clean.append(file_name)
  # get cif file
  cif_file = get_file_from_rcsb(pdb_id,'xray')
  file_to_clean.append(cif_file)
  #
  if os.path.isfile(file_name):
    simple,sym,solv,all = get_cctbx_clashscores(file_name,cif_file=cif_file)
    clashscore_probe = get_probe_clashscore(file_name)
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
    round(sym,1),	            # CCTBX clashscore_only_sym_op
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






