import mmtbx.validation.clashscore
import mmtbx.monomer_library.server
import cctbx.geometry_restraints.manager
from libtbx.utils import Sorry
from libtbx.utils import null_out
from libtbx.utils import Usage
from libtbx import easy_run
from iotbx import pdb
import iotbx.utils
import iotbx.phil
import os,sys
import cProfile
import time



def get_clashscore_probe(file_name, out=sys.stdout):
  '''(str) -> float

  Calculate clashscore using probe.

  Argument:
  file_name: a string, a file name, including the path, of a pdb file
  to which hydorgen were added

  Output:
  clashscore: a float number representning the clashscore of the pdb file
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
    out=out,
    verbose=False)
  return result_with_hydrogen.get_clashscore()

def get_master_phil() :
  # phil parameters for the probe clashscore
  return iotbx.phil.parse("""
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
""")

usage_string = """\
phenix.clashscore file.pdb [params.eff] [options ...]

Options:

  model=input_file      input PDB file
  keep_hydrogens=True   keep input hydrogen files (otherwise regenerate)
  nuclear=False         use nuclear x-H distances and vdW radii
  verbose=True          verbose text output

Example:

  phenix.clashscore model=1ubq.pdb keep_hydrogens=True

"""

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

def get_new_file_name(file_name):
  '''(str) -> str
  When using ready_set, a new file will be created, one that includes
  hydrogens. This function finds the new name and returns it

  Argument:
  file_name: a file name of a pdb file, including the path, that was processed
  using ready_set

  Return:
  file_name: the file name of the modified pdb file, the one that includes hydrogens
  '''
  tmp = file_name.split('.')
  if not (len(tmp)==3 and tmp[1] == 'updated'):
    file_name = get_pdb_file(file_name, print_out=False)
    cmd  = 'phenix.ready_set {}'.format(file_name)
    probe_out = easy_run.go(cmd)
    # getting to modified file name
    line  = [x for x in probe_out.stdout_lines if 'Writing model to' in x]
    if line != []:
      file_name = line[0].split('Writing model to ')[-1]
    else:
      collect_sorry = [x for x in probe_out.stdout_lines if 'Sorry:' in x]
      for x in collect_sorry:
        print x
  return file_name

def get_pdb_file(file_name, print_out=True):
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
  elif print_out:
    print 'Using the file {} found in the working directory'.format(file_name)
  return file_name


if (__name__ == "__main__"):
  file_name = sys.argv[1]
  # get the file name of the file that includes hydrogens
  #file_name = get_new_file_name(file_name)
  #clashscore_probe = get_clashscore_probe(file_name, out=null_out())
  clashscore_probe = get_clashscore_probe(file_name)
  print clashscore_probe
  print 'Done'