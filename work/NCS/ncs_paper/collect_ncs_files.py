"""
Collection of all pdb files with NCS relations, either with only master and
MTRIX records or when NCS is found.

Those files are then filtered for resolution and data requirements as
described in the paper

NCS search if done using the default parameters
"""
from __future__ import division
import sys
import os

__author__ = 'Youval Dar'


class File_records(object):
  """ Information collected on every PDB structure  """

  def __init__(self):

    self.pdb_code = ''
    self.n_ncs_copies = None