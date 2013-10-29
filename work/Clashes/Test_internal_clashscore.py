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
[26,27]
'''

def clashscore(file_name):
  # get the default wroking parameters
  master_phil = get_master_phil()
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
    #log=temp_log)
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
    nonboned_proxies = processed_pdb_file._geometry_restraints_manager._clash_proxies
    clashscore = processed_pdb_file._geometry_restraints_manager._clashscore

def print_record(x):
  res1 = x[0][0][5:-1]
  res2 = x[0][1][5:-1]
  vdw = x[4]
  model = x[3]
  r = model-vdw
  #  N   ILE A 125
  [atom1,resName1,chain1,nRes1] = res1.split()
  [atom2,resName2,chain2,nRes2] = res2.split()
  res1 = '{0:>2} {1:3}  {2:3} {3:>4}'.format(chain1,nRes1, resName1,atom1)
  res2 = '{0:>2} {1:3}  {2:3} {3:>4}'.format(chain2,nRes2, resName2,atom2)
  print '{0}{1} :{2:.3f}'.format(res1,res2,r)

def get_master_phil():

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
  file_name = get_pdb_file(file_name)
  cmd  = 'phenix.ready_set {}'.format(file_name)
  probe_out = easy_run.go(cmd)
  # getting to modified file name
  line  = [x for x in probe_out.stdout_lines if 'Writing model to' in x]
  file_name = line[0].split('Writing model to ')[-1]
  return file_name

def run2 (args, out=sys.stdout) :
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
  set_working_path('Clashes\junk')
  file_name = sys.argv[1]
  # get the file name of the file that includes hydrogens
  file_name = get_new_file_name(file_name)


  run([file_name])