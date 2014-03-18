from __future__ import division
from cctbx.array_family import flex
import os
import mmtbx.model
import libtbx.load_env
from libtbx import easy_pickle
from mmtbx import monomer_library
import mmtbx.monomer_library.server
import mmtbx.monomer_library.pdb_interpretation
from cStringIO import StringIO
from mmtbx import utils
from libtbx.utils import format_cpu_times, null_out
from libtbx.test_utils import approx_equal
import math

def exercise():
  mon_lib_srv = monomer_library.server.server()
  ener_lib = monomer_library.server.ener_lib()
  pdb_file = libtbx.env.find_in_repositories(
                   relative_path="phenix_regression/pdb/enk.pdb", test=os.path.isfile)
  processed_pdb_file = monomer_library.pdb_interpretation.process(
                                       mon_lib_srv               = mon_lib_srv,
                                       ener_lib                  = ener_lib,
                                       file_name                 = pdb_file,
                                       raw_records               = None,
                                       force_symmetry            = True)
  geometry = processed_pdb_file.geometry_restraints_manager(
                                                    show_energies      = False,
                                                    plain_pairs_radius = 5.0)
  #restraints_manager = mmtbx.restraints.manager(geometry      = geometry,
                                                #normalization = False)
  #print '-'*20
  #print dir(geometry) 
  #print '-'*20
  
  #print geometry.pair_proxies
  
  xrs = processed_pdb_file.xray_structure()

  bond_proxies_simple = geometry.pair_proxies(sites_cart =
    xrs.sites_cart()).bond_proxies.simple
  
  sites_cart = xrs.sites_cart()
  
 
  bond_deltas = flex.double()
  rmsdz = flex.double()
  for proxy in bond_proxies_simple:
    #print proxy.i_seqs, proxy.distance_ideal, proxy.weight
    i,j = proxy.i_seqs
    site_1 = sites_cart[i]
    site_2 = sites_cart[j]
    site_deltas = [(x-y)**2 for x,y in zip(site_1,site_2)]
    dist_model = math.sqrt(sum(site_deltas))
    site_bond_deltas = proxy.distance_ideal-dist_model
    sigma = math.sqrt(1.0/proxy.weight)
    z_site = (dist_model - proxy.distance_ideal)/sigma
    bond_deltas.append(site_bond_deltas)
    rmsdz.append(z_site)
    print 'delta bond = %10.4f   z score = %10.4f   sigma = %10.4f'% (site_bond_deltas,z_site,sigma)

    
  STOP()                                              
  xray_structure = processed_pdb_file.xray_structure()
  selection = flex.bool(xray_structure.scatterers().size(), True)
  mol = mmtbx.model.manager(
    restraints_manager = restraints_manager,
    xray_structure = xray_structure,
    pdb_hierarchy = processed_pdb_file.all_chain_proxies.pdb_hierarchy)
  mol.xray_structure.scattering_type_registry(table = "wk1995")


def run():
  exercise()

if (__name__ == "__main__"):
  run()
