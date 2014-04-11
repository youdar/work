from __future__ import division
from iotbx.pdb.multimer_reconstruction import multimer
import mmtbx.refinement.minimization_ncs_constraints
import mmtbx.monomer_library.pdb_interpretation
from scitbx.array_family import flex
import mmtbx.monomer_library.server
from libtbx import adopt_init_args
from libtbx.utils import null_out
from mmtbx import monomer_library
from iotbx.pdb import fetch
import cctbx.adp_restraints
from cctbx import xray
import mmtbx.f_model
import mmtbx.utils
import iotbx.cif
import iotbx.pdb
import os
import sys

"""
Crash when building fmodel for 2qvz.
When removing the SCALE from pdb file it works fine

The following pdb strings were cut from 2qvz
"""


str1 = '''\
CRYST1  631.449  464.724  584.572  90.00 123.84  90.00 C 1 2 1     192
ORIGX1      1.000000  0.000000  0.000000        0.00000
ORIGX2      0.000000  1.000000  0.000000        0.00000
ORIGX3      0.000000  0.000000  1.000000        0.00000
SCALE1      0.001584  0.000000  0.001062        0.00000
SCALE2      0.000000  0.002152  0.000000        0.00000
SCALE3      0.000000  0.000000  0.002060        0.00000
MTRIX1   1  1.000000  0.000000  0.000000        0.00000    1
MTRIX2   1  0.000000  1.000000  0.000000        0.00000    1
MTRIX3   1  0.000000  0.000000  1.000000        0.00000    1
ATOM   3698  N   LEU A 779      97.633  -8.190-322.923  1.00 50.00           N
ATOM   3699  CA  LEU A 779      96.174  -8.136-322.870  1.00 50.00           C
ATOM   3700  C   LEU A 779      95.648  -7.069-323.835  1.00 50.00           C
ATOM   3701  O   LEU A 779      94.524  -6.567-323.614  1.00 50.00           O
ATOM   3702  CB  LEU A 779      95.704  -7.822-321.441  1.00 50.00           C
END
'''

str2 = '''\
CRYST1  631.449  464.724  584.572  90.00 123.84  90.00 C 1 2 1     192
ORIGX1      1.000000  0.000000  0.000000        0.00000
ORIGX2      0.000000  1.000000  0.000000        0.00000
ORIGX3      0.000000  0.000000  1.000000        0.00000
SCALE1      0.001584  0.000000  0.001062        0.00000
SCALE2      0.000000  0.002152  0.000000        0.00000
SCALE3      0.000000  0.000000  0.002060        0.00000
MTRIX1   1  1.000000  0.000000  0.000000        0.00000    1
MTRIX2   1  0.000000  1.000000  0.000000        0.00000    1
MTRIX3   1  0.000000  0.000000  1.000000        0.00000    1
ATOM   3698  N   LEU A 779      97.633  -8.190-322.923  1.00 50.00           N
ATOM   3699  CA  LEU A 779      96.174  -8.136-322.870  1.00 50.00           C
ATOM   3700  C   LEU A 779      95.648  -7.069-323.835  1.00 50.00           C
ATOM   3701  O   LEU A 779      94.524  -6.567-323.614  1.00 50.00           O
ATOM   3702  CB  LEU A 779      95.704  -7.822-321.441  1.00 50.00           C
END
'''

def get_structure_factors():
    """
    Get f_obs and r_free_flags
    From cif file if available
    """
    f_obs = None
    full_path_cif = '2qvz-sf.cif'
    iotbx.cif.reader()
    miller_arrays = iotbx.cif.reader(
      file_path=full_path_cif).\
      as_miller_arrays(force_symmetry=True)
    # print miller_arrays[0].completeness()
    for ma in miller_arrays:
      if ma.is_xray_amplitude_array():
        # Consider using Bijvoet mates
        ma = ma.average_bijvoet_mates()
        f_obs = abs(ma)
        break
      elif not f_obs and ma.is_xray_intensity_array():
        # Consider using Bijvoet mates
        ma = ma.average_bijvoet_mates()
        # convert i_obs to f_obs
        f_obs = abs(ma.french_wilson(log=null_out()))

    if f_obs:
      r_free_flags = f_obs.generate_r_free_flags()
      # f_obs.show_summary()
    else:
      raise RuntimeError("Missing amplitude array.")
    return f_obs,r_free_flags

def run():
  build_fmodel(str1)
  build_fmodel(str2)
  print 'Done...'

def build_fmodel(pdb_str):
  pdb_inp = iotbx.pdb.input(source_info=None, lines=pdb_str)
  ph_ncs = pdb_inp.construct_hierarchy()
  xrs_ncs = pdb_inp.xray_structure_simple()
  f_obs = abs(
    xrs_ncs.structure_factors(d_min=2.0, algorithm="direct").f_calc())
  r_free_flags = f_obs.generate_r_free_flags()

  params = mmtbx.f_model.sf_and_grads_accuracy_master_params.extract()
  params.algorithm = 'direct'
  fmodel = mmtbx.f_model.manager(
    f_obs                        = f_obs,
    r_free_flags                 = r_free_flags,
    xray_structure               = xrs_ncs,
    sf_and_grads_accuracy_params = params,
    target_name                  = 'ml')
  print 'done building fmodel'

if __name__=='__main__':
  os.chdir(r'C:\Phenix\Dev\Work\work\NCS\junk\pdb_test')
  run()


