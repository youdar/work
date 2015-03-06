from __future__ import division
from  iotbx.pdb.multimer_reconstruction import multimer
from libtbx.utils import null_out
from iotbx.pdb import fetch
import iotbx.ncs
import iotbx.pdb
import shutil
import sys
import os

"""
Collection of all pdb files with NCS relations, either with only master and
MTRIX records or when NCS is found.

Those files are then filtered for resolution and data requirements as
described in the paper

NCS search if done using the default NCS search parameters
(minimum chains in master copy, not minimum transforms)
"""


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
    self.only_master_in_pdb = None
    self.ncs_reported_in_pdb = None
    # data_to_param_ration: f_obs.size / atom_number
    self.data_to_param_ration = None
    # refinement_records contains Refinement_results objects
    # for example refinement_records['using cartesian NCS']
    self.refinement_records = {}
    self.r_free_header = None
    self.r_work_header = None
    # list containing issues in data or pdb
    self.issues = []

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
    self.c_beta_deviation = None

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

class ncs_paper_data_collection(object):

  def __init__(self):
    self.files_list = []
    self.ncs_dir = '/net/cci-filer2/raid1/home/youval/work/work/NCS'
    self.asu_dir = self.ncs_dir + '/ncs_paper/ncs_paper_data_files/asu'
    self.mtz_dir = self.ncs_dir + '/ncs_paper/ncs_paper_data_files/mtz'
    self.pdb_dir = self.ncs_dir + '/ncs_paper/ncs_paper_data_files/pdb'
    self.current_dir = os.getcwd()

  def collect_pdb_file(self,pdb_id=None):
    """
    Collect all pdb files with NCS relations and reslution >= 3

    Args:
      pdb_id (str): collect info for a single pdb file
    """
    if not pdb_id:
      # Run on all PDB - Only on LBL machine
      osType = sys.platform
      msg = 'Please run this only on LBL computer'
      assert not osType.startswith('win'),msg
      # set environment
      pdb_dir = os.environ["PDB_MIRROR_PDB"]
      pdb_files = open(os.path.join(pdb_dir, "INDEX"), "r").read().splitlines()
      pdb_codes = [get_4_letters_pdb_id(x) for x in pdb_files]
    else:
      pdb_codes = [pdb_id]

    # get files and save them if have NCS and correct resolution
    for pdb_code in pdb_codes:
      new_rec = File_records()
      new_rec.pdb_code = pdb_code
      fetched_file = pdb_code + '.pdb'
      if not os.path.isfile(fetched_file):
        fetched_file = fetch.fetch(id=pdb_code)
      fn = os.path.realpath(fetched_file)

      pdb_inp = iotbx.pdb.input(file_name=fn)
      mtrix_info = pdb_inp.process_mtrix_records(eps=0.01)
      if len(mtrix_info.r) > 0:
        m = multimer(file_name=fn,reconstruction_type='cau')
        m.write(
          pdb_output_file_name=fn,
          crystal_symmetry=pdb_inp.crystal_symmetry())
        new_rec.n_ncs_copies = m.number_of_transforms + 1
      # note that now the fn is a complete ASU
      ncs_obj = iotbx.ncs.input(file_name=fn)
      if (ncs_obj.number_of_ncs_groups == 0):
        os.remove(fn)
        continue
      new_rec.n_ncs_copies = len(ncs_obj.ncs_transform)
      pio = pdb_inp.get_r_rfree_sigma()
      new_rec.resolution = pio.resolution
      if new_rec.resolution and (new_rec.resolution < 3):
        os.remove(fn)
        continue
      new_rec.year = pdb_inp.extract_header_year()
      new_rec.r_work_header = pio.r_work
      new_rec.r_free_header = pio.r_free
      self.files_list.append(new_rec)
      try:
        shutil.move(fn,self.asu_dir)
      except:
        # avoid error if file already exist
        pass


def get_4_letters_pdb_id(file_name):
  """(str)  -> str
  clean a pdb file name, remove path and file extensions

  Args:
   file_name (str): pdb file name that may look like pdb1a37.pdb
  Return:
   pdb_id (str): the 4 letter pdb id

  >>>get_4_letters_pdb_id('pdb1a37.pdb')
  1a37
  >>>get_4_letters_pdb_id('1a37')
  1a37
  """
  basename = os.path.basename(file_name)
  file_name, file_type = os.path.splitext(basename)
  if len(file_name)>4:
    if 'pdb' in file_name:
      i = file_name.find('pdb')
      pdb_id = file_name[i+3:i+7]
  elif len(file_name)==4:
    pdb_id = file_name
  else:
    pdb_id = None
  return pdb_id
