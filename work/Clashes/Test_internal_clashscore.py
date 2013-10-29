from __future__ import division
from work_func_collection import *
import mmtbx.monomer_library.pdb_interpretation
from mmtbx import monomer_library
import mmtbx.monomer_library.server
import cctbx.geometry_restraints.manager
from libtbx.utils import Sorry
from libtbx import easy_run
import iotbx.utils
import iotbx.phil
import os


''' (str) -> list
Compare clash scores of phenix.clashscore and the internal clashscore function

Argument:
file_name: a pdb file name

Output:
[clashscore_internal,clashscore_probe] : clashscore_internal and clashscore_probe float numbers

>>> python Test_internal_clashscore.py 101m.pdb
101m::26.32::27.10
'''

def clashscore_internal(file_name):
  '''(str) -> float

  Calculate clashscore using pdb_interpertation code.
  
  Argument:
  file_name: a string, a file name, including the path, of a pdb file
  to which hydorgen were added
  
  Output:
  clashscore: a float number representning the clashscore of the pdb file
    
  '''
  # get the default wroking parameters
  master_phil = get_master_phil_interpertation()
  input_objects = iotbx.utils.process_command_line_inputs(
    args=file_name,
    master_phil=master_phil,
    input_types=("pdb","cif"))
  work_phil = master_phil.fetch(sources=input_objects["phil"])
  work_params = work_phil.extract()
  # Get monomers and their properties information
  mon_lib_srv = mmtbx.monomer_library.server.server()
  ener_lib = mmtbx.monomer_library.server.ener_lib()
  # Process pdb file
  processed_pdb_file = mmtbx.monomer_library.pdb_interpretation.process(
    mon_lib_srv=mon_lib_srv,
    ener_lib=ener_lib,
    params=work_params.pdb_interpretation,
    file_name=file_name[0],
    atom_selection_string=work_params.atom_selection,
    strict_conflict_handling=work_params.strict_processing,
    substitute_non_crystallographic_unit_cell_if_necessary=True,
    max_atoms=work_params.max_atoms,
    log=sys.stdout)
  # conflict handling
  if (work_params.strict_processing):
    msg = processed_pdb_file.all_chain_proxies.fatal_problems_message()
    if (msg is not None):
      raise Sorry(msg)


  if (work_params.build_geometry_restraints_manager):
    processed_pdb_file.geometry_restraints_manager(
      params_edits=work_params.geometry_restraints.edits,
      params_remove=work_params.geometry_restraints.remove)
    site_cart = processed_pdb_file._geometry_restraints_manager.sites_cart_used_for_pair_proxies()
    # Get clashscore
    return processed_pdb_file._geometry_restraints_manager._clashscore
  else:
    # a clashscore of -1 indicates the process did not work
    return -1


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
  file_name = get_pdb_file(file_name)
  cmd  = 'phenix.ready_set {}'.format(file_name)
  probe_out = easy_run.go(cmd)
  # getting to modified file name
  line  = [x for x in probe_out.stdout_lines if 'Writing model to' in x]
  file_name = line[0].split('Writing model to ')[-1]
  return file_name

def clashscore_probe(file_name, out=sys.stdout):
  '''(str) -> float

  Calculate clashscore using probe.
  
  Argument:
  file_name: a string, a file name, including the path, of a pdb file
  to which hydorgen were added
  
  Output:
  clashscore: a float number representning the clashscore of the pdb file
  '''
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


if (__name__ == "__main__"):
  set_working_path('Clashes\junk')
  file_name = sys.argv[1]
  # get the file name of the file that includes hydrogens
  file_name = get_new_file_name(file_name)
  #
  print 'Starting'
  print os.getcwd()
  print file_name
  print '*'*80

  clashscore_internal = clashscore_internal(file_name)
  clashscore_probe = clashscore_probe(file_name, out=sys.stdout)
  
  print '*'*80
  print 'done'
