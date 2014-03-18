from __future__ import division
from iotbx.pdb.multimer_reconstruction import multimer
import mmtbx.refinement.minimization_ncs_constraints
import mmtbx.monomer_library.pdb_interpretation
from libtbx.test_utils import approx_equal
from scitbx.array_family import flex
import mmtbx.monomer_library.server
from libtbx import adopt_init_args
from mmtbx import monomer_library
from cctbx import xray
import mmtbx.f_model
import mmtbx.utils
import iotbx.pdb
import os
import sys

class ncs_refine_test(object):
  def __init__(self,
               n_macro_cycle,
               sites,
               u_iso,
               finite_grad_differences_test,
               use_geometry_restraints,
               d_min):
    """
    Arguments:
    n_macro_cycle:  number of refinement cycles


    """

  def refine_using_strict_ncs(self):
    """

    """

  def refine_using_complete_asu(self):
    """

    """