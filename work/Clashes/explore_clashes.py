from __future__ import division
from work_func_collection import *
import mmtbx.monomer_library.pdb_interpretation
from mmtbx import monomer_library
import mmtbx.monomer_library.server
import cctbx.geometry_restraints.manager
from libtbx.utils import Sorry
from libtbx.utils import null_out
from libtbx import easy_run
import iotbx.utils
import iotbx.phil
import os


'''
This code a modification of code taken from
c:\Phenix\Dev\phenix_sources\cctbx_project\mmtbx\command_line\pdb_interpretation.py
'''

def run(file_name):
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
    #log=sys.stdout)
    log=null_out())
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
    # get_sorted(self, by_value, sites_cart, site_labels=None, max_items=None)
    #nonboned_proxies = processed_pdb_file._geometry_restraints_manager._pair_proxies.nonbonded_proxies.get_sorted(
      #by_value="delta",
      #sites_cart=site_cart)
    nonboned_proxies = processed_pdb_file._geometry_restraints_manager._clash_proxies
    clashscore = processed_pdb_file._geometry_restraints_manager._clashscore
    print '*'*80
    print '*'*80
    for i,x in enumerate(nonboned_proxies):
      print_record(x)
      if i>19: break
    print '*'*80
    print 'Clash list length: {}'.format(len(nonboned_proxies))
    print 'clashscore: {0:.2f}'.format(clashscore)
    print '*'*80
    #print len(nonboned_proxies[0])
    #print clashscore
    #print len(nonboned_proxies[0])/clashscore

    #deltas = [x[3] for x in nonboned_proxies[0]]
    #print 'Min delta: {0:.3f}'.format(min(deltas))
    #print 'Max delta: {0:.3f}'.format(max(deltas))
    #print '*'*70
    #i = 0
    #for x in nonboned_proxies[0]:
      #res1 = x[0][0][5:-1]
      #res2 = x[0][1][5:-1]
      #vdw = x[4]
      #model = x[3]
      #r = model-vdw
      #if r<=-0.4:
        #i += 1
        #print '{0} {1}  model= {2:.3f} vdw = {3:.3f} model-vdwl= {4:.3f}'.format(res1,res2,model,vdw,r)
      #if ('30' in res1) and ('30' in res2):
        #print '*'*70
        #print '{0} {1}  model= {2:.3f} vdw = {3:.3f} model-vdwl= {4:.3f}'.format(res1,res2,model,vdw,r)
    #print '*'*70
    #print 'number of clashes: {}'.format(i)

    #for i in range(10):
      #res1 = nonboned_proxies[0][i][0][0][5:-1]
      #res2 = nonboned_proxies[0][i][0][1][5:-1]
      #delta = nonboned_proxies[0][i][3]
      #print '{0} {1}  :  {2:.3f}'.format(res1,res2,delta)

    print 'ok'

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



if __name__=='__main__':
  print '*'*80
  # set work folder
  set_working_path('Clashes\junk')
  print 'Current working path {}'.format(os.getcwd())
  #file_name = get_pdb_file('3e8k.pdb')
  file_name = get_pdb_file('101m_H.pdb')
  #cmd  = 'phenix.ready_set {}'.format(file_name)
  #print 'Running: {}'.format(cmd)
  #print '*'*80
  #probe_out = easy_run.go(cmd)
  ## getting to modified file name
  #line  = [x for x in probe_out.stdout_lines if 'Writing model to' in x]
  #file_name = line[0].split('Writing model to ')[-1]
  print 'Using the file {}, the file returned from phenix.ready_set'.format(file_name)
  #print '*'*80
  out = easy_run.go('phenix.clashscore {}  keep_hydrogens=True'.format(file_name))
  #print 'phenix.{}'.format(out.stdout_lines[-1])
  for l in out.stdout_lines:
    print l
  print '*'*80
  # Call the new clashscore
  run([file_name])