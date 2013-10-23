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
    log=sys.stdout)
    #log=None)
  # conflict handling 
  if (work_params.strict_processing):
    msg = processed_pdb_file.all_chain_proxies.fatal_problems_message()
    if (msg is not None):
      raise Sorry(msg)
    
  # processed_pdb_file.geometry_restraints_manager options
  # geometry_restraints_manager(self, plain_pairs_radius=None, params_edits=None, params_remove=None, 
  # hydrogen_bond_proxies=None, hydrogen_bond_params=None, custom_nonbonded_exclusions=None, 
  # assume_hydrogens_all_missing=True, show_energies=True, hard_minimum_bond_distance_model=0.001, 
  # external_energy_function=None, ramachandran_atom_selection=None, den_manager=None, 
  # reference_manager=None) 

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

    print len(nonboned_proxies[0])
    print clashscore
    print len(nonboned_proxies[0])/clashscore
    
    deltas = [x[3] for x in nonboned_proxies[0]]
    print 'Min delta: {0:.3f}'.format(min(deltas))
    print 'Max delta: {0:.3f}'.format(max(deltas))
    print '*'*70
    i = 0
    for x in nonboned_proxies[0]:
      res1 = x[0][0][5:-1]
      res2 = x[0][1][5:-1]
      vdw = x[4]
      model = x[3]
      r = model-vdw 
      if r<=-0.4:
        i += 1
        print '{0} {1}  model= {2:.3f} vdw = {3:.3f} model-vdwl= {4:.3f}'.format(res1,res2,model,vdw,r)
      if ('30' in res1) and ('30' in res2):
        print '*'*70
        print '{0} {1}  model= {2:.3f} vdw = {3:.3f} model-vdwl= {4:.3f}'.format(res1,res2,model,vdw,r)
    print '*'*70      
    print 'number of clashes: {}'.format(i)
     
    #for i in range(10):  
      #res1 = nonboned_proxies[0][i][0][0][5:-1]
      #res2 = nonboned_proxies[0][i][0][1][5:-1]
      #delta = nonboned_proxies[0][i][3]
      #print '{0} {1}  :  {2:.3f}'.format(res1,res2,delta)
    
    print 'ok'

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
  # set work folder
  set_working_path('Clashes\junk')
  print 'Current working path {}'.format(os.getcwd())
  file_name = [get_pdb_file('3e8k')]
  run(file_name)