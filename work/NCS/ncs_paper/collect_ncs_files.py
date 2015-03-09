from __future__ import division
from  iotbx.pdb.multimer_reconstruction import multimer
from libtbx.utils import null_out
from iotbx.pdb import fetch
import cPickle as pickle
import iotbx.ncs
import iotbx.pdb
import shutil
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

    self.pdb_id = ''
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
    self.map_cc = None
    self.rmsd = None
    self.rama_outliers = None
    self.rotamer_outliers = None

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

  def get_pdb_file_info(self,pdb_id):
    """
    Collect pdb file NCS and header info if there are NCS  relations
    and reslution >= 3

    Args:
      pdb_id (str): collect info for a single pdb file

    Return:
      new_rec (File_records object): object containing the collected info
    """
    new_rec = File_records()
    new_rec.pdb_id = pdb_id
    fetched_file = pdb_id + '.pdb'
    if not os.path.isfile(fetched_file):
      fetched_file = fetch.get_pdb (
        id=pdb_id, data_type='pdb',
        mirror='rcsb',quiet=True,log=null_out())
      # pdb_obj = fetch.fetch(id=pdb_id)
      # pdb_id = get_4_letters_pdb_id(pdb_obj.filename)
      # fetched_file = pdb_id + '.pdb'
      # pdb_obj.write(fetched_file)
    else:
      print 'file exist in local folder, did not fetch it..'
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
      return None
    new_rec.n_ncs_copies = len(ncs_obj.ncs_transform)
    pio = pdb_inp.get_r_rfree_sigma()
    new_rec.resolution = pio.resolution
    if new_rec.resolution and (new_rec.resolution < 3):
      os.remove(fn)
      return None
    new_rec.year = pdb_inp.extract_header_year()
    new_rec.r_work_header = pio.r_work
    new_rec.r_free_header = pio.r_free
    try:
      shutil.move(fn,self.asu_dir)
    except:
      # avoid error if file already exist
      pass
    return new_rec

  def write_to_file(self,file_name,file_record):
    """
    writes pickled object to file_name
    """
    if file_record:
      pickle.dump(file_record,open(file_name,'w'))

  def make_mtz_file(self):
    """ get cif file and create mtz """
    pass

def get_cif_file(self,pdb_id, write_folder):
  """ get the cif file for pdb_id and write it in write_folder"""
  fetched_file = fetch.get_pdb (
    id=self.pdb_code, data_type='xray',
    mirror='rcsb',quiet=True,log=null_out())
  if os.path.isdir(write_folder):
    try:
      shutil.move(fetched_file,write_folder)
      return True
    except:
      # avoid error if file already exist
      pass
  return False


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

def get_structure_factors(self):
    """
    Get f_obs and r_free_flags From cif file, if available
    """
    try:
      inputs = mmtbx.utils.process_command_line_args(args = self.args)
      df = mmtbx.utils.determine_data_and_flags(
        reflection_file_server  = inputs.get_reflection_file_server(),
        log = null_out())
      self.f_obs = df.f_obs
      self.r_free_flags = df.r_free_flags
    except:
      if self.print_during_refinement:
        print 'Getting f_obs and r_free_flags by processing cif file'
      self.f_obs = None
      self.r_free_flags = None
      fobs = ["FOBS,SIGFOBS",'FOBS','FOBS,PHIM',
              "F(+),SIGF(+),F(-),SIGF(-)","F(+),F(-)"]
      iobs = ["IOBS,SIGIOBS",'IOBS','IOBS,PHIM',
              'I(+),SIGI(+),I(-),SIGI(-)']
      self.i_obs = None
      self.f_obs = None
      if self.full_path_cif:
        miller_arrays = self.get_miller_arrays()
        # print miller_arrays[0].completeness()
        for ma in miller_arrays:
          ls = ma.info().label_string()
          if (ls in fobs):
            # Consider using Bijvoet mates
            ma = ma.average_bijvoet_mates()
            self.f_obs = abs(ma)
          elif ls == "R-free-flags":
            self.r_free_flags = abs(ma)
          elif not self.f_obs and (ls in iobs):
            # Consider using Bijvoet mates
            ma = ma.average_bijvoet_mates()
            # convert i_obs to f_obs
            self.i_obs = ma
          elif not self.r_free_flags and ls == "R-free-flags-1":
            self.r_free_flags = abs(ma.french_wilson(log=null_out()))
      else:
        raise RuntimeError("No cif file.")

      # When fobs where not found via string look of other fobs forms
      if not self.f_obs:
        for ma in miller_arrays:
          if ma.is_xray_amplitude_array():
            ma = ma.average_bijvoet_mates()
            self.f_obs = abs(ma)
      if (not self.f_obs) and (not self.i_obs):
        for ma in miller_arrays:
          if ma.is_xray_intensity_array():
            ma = ma.average_bijvoet_mates()
            self.i_obs = ma
      if not self.f_obs and self.i_obs:
        self.f_obs = abs(self.i_obs.french_wilson(log=null_out()))
    #
    if self.f_obs:
      # self.f_obs.show_summary()
      self.data_set_size = self.f_obs.size()
      if self.r_free_flags:
        self.f_obs, self.r_free_flags = self.f_obs.common_sets(
          self.r_free_flags)
        self.r_free_flags = self.make_r_free_boolean(self.r_free_flags)
      else:
        self.r_free_flags = self.f_obs.generate_r_free_flags()
      # Data completeness: Fraction of unmeasured reflections within the
      # [d_min, d_max] range,where d_min and d_max are highest and lowest
      # resolution of data set correspondingly.
      self.completeness = self.f_obs.array().completeness()
    else:
      raise RuntimeError("Missing amplitude array.")
