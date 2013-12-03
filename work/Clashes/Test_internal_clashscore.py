from __future__ import division
import mmtbx.monomer_library.pdb_interpretation
from mmtbx import monomer_library
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

''' (str) -> str
Compare clash scores of phenix.clashscore and the internal clashscore function

Argument:
file_name: a pdb file name

Output:
'pdb_file_name::clashscore_internal::clashscore_probe'

>>> python Test_internal_clashscore.py 3e8k.pdb
3e8k::35.210::62.472
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
  pdb = monomer_library.pdb_interpretation.run(args=[file_name],
    substitute_non_crystallographic_unit_cell_if_necessary=True,
    assume_hydrogens_all_missing=False,
    hard_minimum_nonbonded_distance=0.0,
    log=null_out())

  grm = pdb.geometry_restraints_manager(assume_hydrogens_all_missing=False,
                                        hard_minimum_nonbonded_distance=0.0)

  xrs = pdb.xray_structure()
  sites_cart = xrs.sites_cart()
  site_labels = pdb._geometry_restraints_manager._site_lables
  hd_sel = xrs.hd_selection()
  table_bonds = grm.shell_sym_tables[0]
  full_connectivty_table = table_bonds.full_simple_connectivity()

  grm.get_nonbonded_clashscore(sites_cart=sites_cart,
                               site_labels=site_labels,
                               hd_sel=hd_sel,
                               full_connectivty_table=full_connectivty_table)

  clashscore_all = grm.nonbonded_clash_info.nb_clashscore_all_clashes
  clashscore_without_sym_op = grm.nonbonded_clash_info.nb_clashscore_without_sym_op
  clashscore_due_to_sym_op = grm.nonbonded_clash_info.nb_clashscore_due_to_sym_op

  return clashscore_all,clashscore_without_sym_op,clashscore_due_to_sym_op

  ## get the default wroking parameters
  #master_phil = get_master_phil_interpertation()
  #input_objects = iotbx.utils.process_command_line_inputs(
    #args=[file_name],
    #master_phil=master_phil,
    #input_types=("pdb","cif"))
  #work_phil = master_phil.fetch(sources=input_objects["phil"])
  #work_params = work_phil.extract()
  ## Get monomers and their properties information
  #mon_lib_srv = mmtbx.monomer_library.server.server()
  ## Get bonds properties from
  ## phenix_sources\cctbx_project\mmtbx\monomer_library\server.py
  ## phenix_sources/chem_data/mon_lib/list/mon_lib_list.cif

  #ener_lib = mmtbx.monomer_library.server.ener_lib()
  ## Process pdb file
  #processed_pdb_file = mmtbx.monomer_library.pdb_interpretation.process(
    #mon_lib_srv=mon_lib_srv,
    #ener_lib=ener_lib,
    #params=work_params.pdb_interpretation,
    #file_name=file_name,
    #atom_selection_string=work_params.atom_selection,
    #strict_conflict_handling=work_params.strict_processing,
    #substitute_non_crystallographic_unit_cell_if_necessary=True,
    #max_atoms=work_params.max_atoms,
    #log=null_out())
    ##log=sys.stdout)
  ## conflict handling
  #if (work_params.strict_processing):
    #msg = processed_pdb_file.all_chain_proxies.fatal_problems_message()
    #if (msg is not None):
      #raise Sorry(msg)


  #if (work_params.build_geometry_restraints_manager):
    #processed_pdb_file.geometry_restraints_manager(
      #assume_hydrogens_all_missing=False,
      #hard_minimum_nonbonded_distance=0.0,
      #params_edits=work_params.geometry_restraints.edits,
      #params_remove=work_params.geometry_restraints.remove)
    ## Get clashscore
    #clashscore_all = processed_pdb_file.nonbonded_clash_info.nb_clashscore_all_clashes
    #clashscore_only_sym_op = processed_pdb_file.nonbonded_clash_info.nb_clashscore_only_sym_op
    #clashscore_without_sym_op = processed_pdb_file.nonbonded_clash_info.nb_clashscore_without_sym_op
    #return clashscore_all,clashscore_without_sym_op,clashscore_only_sym_op
  #else:
    ## a clashscore of -1 indicates the process did not work
    #return -1



def tic():
  #Homemade version of matlab tic and toc functions
  global startTime_for_tictoc
  startTime_for_tictoc = time.time()

def toc(msg='',print_time=True):
  if 'startTime_for_tictoc' in globals():
    if print_time:
      outstr = '{0}: Elapsed time is: {1:.4f} seconds\n'.format(msg,time.time() - startTime_for_tictoc)
      print outstr
    else:
      outstr = '{0:.4f}'.format(time.time() - startTime_for_tictoc)
      return outstr
  else:
    print "Toc: start time not set"

def get_master_phil_interpertation():
  # this is a phil for the pdb_interpertation execution

  return iotbx.phil.parse("""\
atom_selection = None
  .type = str
  .help = "Limit all analysis of restraints to this selection only."
strict_processing = False
  .type = bool
build_geometry_restraints_manager = True
  .type = bool
build_xray_structure = True
  .type = bool
max_atoms = None
  .type = int

write_geo_files = False
  .type = bool
write_tardy_geo_files = False
  .type = bool


include scope mmtbx.monomer_library.pdb_interpretation.grand_master_phil_str
""", process_includes=True)



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
  # We to keep hydrogen
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

def call_both_clashscores(file_name):
  tic()
  clashscore_all,clashscore_without_sym_op,clashscore_only_sym_op = get_clashscore_internal(file_name)
  time_internal = toc(print_time=False)
  #clashscore_all = clashscore_without_sym_op = clashscore_only_sym_op = 0
  clashscore_internal = [clashscore_all,clashscore_without_sym_op,clashscore_only_sym_op]
  tic()
  #clashscore_probe = get_clashscore_probe(file_name, out=null_out())
  #clashscore_probe = get_clashscore_probe(file_name)
  time_probe = toc(print_time=False)
  clashscore_probe = 0
  return clashscore_internal,clashscore_probe,time_internal,time_probe


if (__name__ == "__main__"):
  set_working_path('Clashes\wtest')
  file_name = sys.argv[1]
  # get the file name of the file that includes hydrogens
  file_name = get_new_file_name(file_name)
  #
  #print 'Starting'
  #print os.getcwd()
  #print file_name
  #print '*'*80


  # Performance test
  #cProfile.run("clashscore_internal,clashscore_probe = call_both_clashscores(file_name)")
  #print '*'*80
  #print 'done'

  # regular run
  clashscore_internal,clashscore_probe,time_internal,time_probe = call_both_clashscores(file_name)

  #print '\nClashscore all: {0:.2f}\nwithout_sym_op: {1:.2f}\nonly_sym_op   : {2:.3f}\n'.format(*clashscore_internal)
  output_file_name = get_file_name(file_name)
  outstr = '{0}::{1:.1f}::{2:.1f}::{3}::{4}'.format(output_file_name,clashscore_internal[0],clashscore_probe,time_internal,time_probe)
  print outstr




