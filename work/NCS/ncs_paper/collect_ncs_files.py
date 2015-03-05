"""
Collection of all pdb files with NCS relations, either with only master and
MTRIX records or when NCS is found.

Those files are then filtered for resolution and data requirements as
described in the paper

NCS search if done using the default NCS search parameters
(minimum chains in master copy, not minimum transforms)
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
    self.year = None
    self.resolution = None
    self.data_completeness = None
    self.solvent = None
    # data_to_param_ration: f_obs.size / atom_number
    self.data_to_param_ration = None
    # refinement_records contains Refinement_results objects
    self.refinement_records = {}
    self.r_free_header = None
    self.r_work_header = None

  def __repr__(self):
    """ prints object's summary  """
    s = '{:<35}:{:<10}'
    out_lst = get_dict_as_list(self.__dict__,s)
    return '\n'.join(out_lst)

class Refinement_results(object):
  """
  Collect the results of a particular refinement test
  """

  def __init__(self):
    self.r_free_init = 0
    self.r_work_init = 0
    self.r_free_final = 0
    self.r_work_final = 0
    self.refinement_time = None
    self.normalized_sym_nbo = None
    self.clashscore = None

  def __repr__(self):
    """ prints object's summary  """
    s = '{:<35}:{:<10}'
    out_lst = get_dict_as_list(self.__dict__,s)
    out_lst.sort()
    return '\n'.join(out_lst)


def get_dict_as_list(d,template):
  """
  recursively expands dictionary for printing

  Args:
    d (dict):  a dictionary
    template (str): a template used to format the key and value of dict

  Returns:
    out_lst (list): a list of string containing the formatted key-value pairs
  """
  out_lst = []
  keys = d.keys()
  keys.sort()
  for k in keys:
    v = d[k]
    if type(v) is dict:
      x = get_dict_as_list(v,template)
      out_lst.extend(x)
    elif 'Refinement_results' in type(v).__name__:
      out_lst.append(v.__repr__())
    elif v and (v !=  0):
      out_lst.append(template.format(k,v))
  return out_lst

