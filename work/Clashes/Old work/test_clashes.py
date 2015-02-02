from __future__ import division
import mmtbx.validation.clashscore
from mmtbx.tls.tools import tls_from_pdb_inp
from iotbx import pdb
import iotbx.phil
from libtbx.utils import Usage
import cProfile
import sys,os


def run (args, out=sys.stdout) :
  '''(str) -> str

  This function tests the clash scores of pdb files

  It returns pdb_file_name::score_with_hydrogen::score_without_hydrogen::experment_type

  >>> python test_clashes.py 4iw4
  4iw4::8.9174::8.9174::X-RAY DIFFRACTION

  '''
  file_name = args[0]
  if not os.path.isfile(file_name):
    osType = sys.platform
    if osType.startswith('win'):
      file_name = get_pdb_file(file_name)
    else:
      file_name = get_pdb_file_dir(file_name)
    args[0] = file_name
  experment_type = get_experment_type(file_name)
  cmdline = iotbx.phil.process_command_line_with_files(
    args=args,
    master_phil=get_master_phil(),
    pdb_file_def="model",
    usage_string=usage_string)
  params = cmdline.work.extract()
  if (params.model is None) :
    raise Usage(usage_string)
  pdb_in = cmdline.get_file(params.model, force_type="pdb")
  hierarchy = pdb_in.file_object.construct_hierarchy()
  # We want to check results with and without hydrogen
  result_with_hydrogen = mmtbx.validation.clashscore.clashscore(
    pdb_hierarchy=hierarchy,
    keep_hydrogens=True,
    nuclear=params.nuclear,
    out=out,
    verbose=False)
  result_without_hydrogen = mmtbx.validation.clashscore.clashscore(
      pdb_hierarchy=hierarchy,
      keep_hydrogens=False,
      nuclear=params.nuclear,
      out=out,
      verbose=False)
  # build output string
  file_name  = get_file_name(file_name)
  score_with_h = result_with_hydrogen.get_clashscore()
  score_without_h = result_without_hydrogen.get_clashscore()
  print '{0}::{1:.4f}::{2:.4f}::{3}'.format(file_name,score_with_h,score_without_h,experment_type)

def get_master_phil () :
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

def get_pdb_file(file_name):
  from iotbx.pdb import fetch
  class null_out(object):
    """Pseudo-filehandle for suppressing printed output."""
    def isatty(self): return False
    def close(self): pass
    def flush(self): pass
    def write(self, str): pass
    def writelines(self, sequence): pass

  log  = null_out()
  file_name = get_file_name(file_name)
  return fetch.get_pdb (file_name,'pdb',mirror='rcsb',log=log)

def get_pdb_file_dir(file_name):
  file_name = get_file_name(file_name)
  pdb_dir = os.environ["PDB_MIRROR_PDB"]
  pdb_files = open(os.path.join(pdb_dir, "INDEX"), "r").read().splitlines()
  for i,p in enumerate(pdb_files):
    if file_name in p:
      break
  return os.path.join(pdb_dir,pdb_files[i])


def set_working_path():
  # locate the directory containing the log files
  osType = sys.platform
  if osType.startswith('win'):
      directory_path = 'c:\Phenix\Dev\Work\work\Clashes'
  else:
      directory_path = '/net/cci-filer2/raid1/home/youval/Work/work/Clashes'
  os.chdir(directory_path)

def get_file_name(file_name):
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

def get_experment_type(file_name):
  '''(str) -> str
  Look for EXPERIMENT TYPE in PDB REMARK
  REMARK 200 		: X-RAY DIFFRACTION
  REMARK 210,215,217 	: NMR
  REMARK 230		: NEUTRON DIFFRACTION
  REMARK 245		: ELECTRON MICROSCOPE
  REMARK 250		: Other
  REMARK 265		: SMALL ANGLE X-RAY SCATTERING

  returns a string combined from all experiment types

  Note that the experiment identification is not done by actually reading the
  records, but by the record clasification at
  http://www.wwpdb.org/documentation/format33/remarks1.html
  '''
  remark_dict = dict([(200,'X-RAY DIFFRACTION'),
                      (210,'NMR'),(215,'NMR'),(217,'NMR'),
                      (230,'NEUTRON DIFFRACTION'),
                      (245,'ELECTRON MICROSCOPE'),
                      (250,'Other'),
                      (265,'SMALL ANGLE X-RAY SCATTERING')])
  experment_type = []
  for iii in remark_dict:
    rem_records = get_remark_iii_records(file_name,iii)
    if rem_records != []:
      experment_type.append(remark_dict[iii])

  return '::'.join(experment_type)

def get_remark_iii_records(file_name,iii):
  pdb_inp = pdb.hierarchy.input(file_name=file_name)
  pdb_inp_tls = tls_from_pdb_inp(
    remark_3_records = pdb_inp.input.extract_remark_iii_records(iii),
    pdb_hierarchy = pdb_inp.hierarchy)
  return pdb_inp_tls.remark_3_records

if (__name__ == "__main__"):
  set_working_path()
  run(sys.argv[1:])
