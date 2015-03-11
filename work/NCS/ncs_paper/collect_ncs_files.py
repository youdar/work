from __future__ import division
from  iotbx.pdb.multimer_reconstruction import multimer
from iotbx.ncs.ncs_preprocess import ncs_only
import mmtbx.monomer_library.server
import iotbx.reflection_file_utils
from libtbx.utils import null_out
from iotbx.pdb import fetch
from libtbx import easy_run
import cPickle as pickle
from glob import glob
import mmtbx.utils
import mmtbx.masks
import iotbx.ncs
import iotbx.pdb
import iotbx.mtz
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

    self.pdb_id = ''
    self.n_ncs_copies = None
    self.year = None
    self.resolution = None
    self.data_completeness = None
    self.solvent_fraction = None
    self.experiment_type = None
    self.only_master_in_pdb = None
    self.ncs_reported_in_pdb = None
    self.n_atoms_in_asu = None
    # data_to_param_ratio: f_obs.size / atom_number
    self.data_to_param_ratio = None
    self.data_to_param_ratio_ncs = None
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

class ncs_paper_data_collection(object):

  def __init__(self):
    self.files_list = []
    self.pdbs_dict = {}
    osType = sys.platform
    if osType.startswith('win'):
      s = r'C:\Phenix\Dev\work\work\work\NCS\ncs_paper\ncs_paper_data_files'
      self.ncs_dir = s
    else:
      s= '/net/cci/youval/work/work/NCS/ncs_paper/ncs_paper_data_files/'
      self.ncs_dir = s
    self.asu_dir = os.path.join(self.ncs_dir,'asu')
    self.mtz_dir = os.path.join(self.ncs_dir,'mtz')
    self.pdb_dir = os.path.join(self.ncs_dir,'pdb')
    self.data_dir = os.path.join(self.ncs_dir,'data')
    self.model_vs_data_dir = os.path.join(self.ncs_dir,'model_vs_data')
    self.pdb_not_used_dir = os.path.join(self.ncs_dir,'pdb_with_ncs_not_used')
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
    else:
      print 'file exist in local folder, did not fetch it..'
    fn = os.path.realpath(fetched_file)
    pdb_inp = iotbx.pdb.input(file_name=fn)
    mtrix_info = pdb_inp.process_mtrix_records(eps=0.01)
    t = (mtrix_info.as_pdb_string() == '')
    t |= (not ncs_only(mtrix_info))
    new_rec.only_master_in_pdb = not t
    if len(mtrix_info.r) > 0 and new_rec.only_master_in_pdb:
      m = multimer(file_name=fn,reconstruction_type='cau')
      m.write(
        pdb_output_file_name=fn,
        crystal_symmetry=pdb_inp.crystal_symmetry())
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
    new_rec.experiment_type = pdb_inp.get_experiment_type()
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

  def read_from_file(self,file_name,path=''):
    """
    reads pickled object from file_name

    Args:
      file_name (str): such as log_"pdb id"

    Returns:
      File_records object
    """
    fn = os.path.join(path,file_name)
    if os.path.isfile(fn):
      return pickle.load(open(fn,'r'))
    else:
      return None

  def make_mtz_file(self,file_record):
    """
    get cif file and create mtz and add info to pdb_file_records

    Args:
      file_record (obj): File_records object

    Return:
      updated file_record
    """
    cif = get_cif_file(file_record.pdb_id)
    pdb = os.path.join(self.asu_dir,file_record.pdb_id + '.pdb')
    if cif:
      f_obs,r_free_flags,completeness,data_size = get_structure_factors(
        pdb,cif,self.mtz_dir)
    else:
      return None
    if f_obs:
      file_record.data_completeness = completeness
      pdb_inp = iotbx.pdb.input(file_name=pdb)
      #
      n_atoms = pdb_inp.atoms().size()
      file_record.data_to_param_ratio = data_size/3.0/n_atoms
      file_record.n_atoms_in_asu = n_atoms
      #
      xrs_asu = pdb_inp.xray_structure_simple()
      file_record.solvent_fraction = mmtbx.masks.asu_mask(
        xray_structure=xrs_asu,
        d_min=f_obs.d_min()).asu_mask.contact_surface_fraction
      return file_record
    else:
      return None

  def collect_all_file_records(self):
    """
    Collect the information on all files

    Returns a dictionary of all data collected
      dict key: pdb_id
      dict val: File_records object
    """
    files_list = glob(self.data_dir + '/log_*')
    for fn in files_list:
      r = pickle.load(open(fn,'r'))
      self.pdbs_dict[r.pdb_id] = r
    self.files_list = [x[-4:] for x in files_list]
    return self.pdbs_dict


def get_cif_file(pdb_id):
  """ get the cif file """
  try:
    fetched_file = fetch.get_pdb(
      id=pdb_id, data_type='xray',
      mirror='rcsb',quiet=True,log=null_out())
    return fetched_file
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
    else:
      pdb_id = file_name.replace('-sf','')
  elif len(file_name)==4:
    pdb_id = file_name
  else:
    pdb_id = None
  return pdb_id

def get_structure_factors(pdb,cif,mtz_folder):
  """
  Get f_obs and r_free_flags From cif file

  Args:
    mtz_folder (str): path to the folder mtz file will be saved to
    pdb (str): pdb file path
    cif (str): cif file path

  Returns:
    f_obs, i_obs, r_free_flags, completeness, data_set_size
      Data completeness: Fraction of unmeasured reflections within the
  """
  f_obs = None
  i_obs = None
  r_free_flags = None
  if not (os.path.isfile(pdb) and os.path.isfile(cif)):
    return None,None,0,0
  miller_arrays = get_miller_arrays(pdb,cif,mtz_folder)
  try:
    inputs = mmtbx.utils.process_command_line_args(args = [pdb,cif])
    df = mmtbx.utils.determine_data_and_flags(
      reflection_file_server  = inputs.get_reflection_file_server(),
      log = null_out())
    f_obs = df.f_obs
    r_free_flags = df.r_free_flags
  except:
    # if the simple way did not work try the following
    fobs_type = ["FOBS,SIGFOBS",'FOBS','FOBS,PHIM',
                 "F(+),SIGF(+),F(-),SIGF(-)","F(+),F(-)"]
    iobs_type = ["IOBS,SIGIOBS",'IOBS','IOBS,PHIM',
                 'I(+),SIGI(+),I(-),SIGI(-)']
    # print miller_arrays[0].completeness()
    if miller_arrays:
      for ma in miller_arrays:
        ls = ma.info().label_string()
        if (ls in fobs_type):
          # Consider using Bijvoet mates
          ma = ma.average_bijvoet_mates()
          f_obs = abs(ma)
        elif ls == "R-free-flags":
          r_free_flags = abs(ma)
        elif not f_obs and (ls in iobs_type):
          # Consider using Bijvoet mates
          ma = ma.average_bijvoet_mates()
          # convert i_obs to f_obs
          i_obs = ma
        elif not r_free_flags and ls == "R-free-flags-1":
          r_free_flags = abs(ma.french_wilson(log=null_out()))

      # When fobs where not found via string look of other fobs forms
      if (not f_obs) and (not i_obs):
        for ma in miller_arrays:
          if ma.is_xray_amplitude_array():
            ma = ma.average_bijvoet_mates()
            f_obs = abs(ma)
      if (not f_obs) and (not i_obs):
        for ma in miller_arrays:
          if ma.is_xray_intensity_array():
            ma = ma.average_bijvoet_mates()
            i_obs = ma
      if not f_obs and i_obs:
        f_obs = abs(i_obs.french_wilson(log=null_out()))
  #
  if f_obs:
    # f_obs.show_summary()
    data_set_size = f_obs.size()
    if r_free_flags:
      f_obs, r_free_flags = f_obs.common_sets(r_free_flags)
      r_free_flags = make_r_free_boolean(r_free_flags)
    else:
      r_free_flags = f_obs.generate_r_free_flags()
    # Data completeness: Fraction of unmeasured reflections within the
    # [d_min, d_max] range,where d_min and d_max are highest and lowest
    # resolution of data set correspondingly.
    completeness = f_obs.array().completeness()
  else:
    return None,None,0,0
  return f_obs, r_free_flags, completeness,data_set_size

def get_miller_arrays(pdb,cif,mtz_folder):
  """
  convert cif to mtz and write it in the mtz_folder

  convert cif of the format 'r_name_sf.ent.gz' to mtz file
  creates mtz file with crystal symmetry in current folder

  Returns:

  """
  if not (os.path.isfile(pdb) and os.path.isfile(cif)):
    return None
  pdb_id = get_4_letters_pdb_id(cif)
  mtz_fn = os.path.join(mtz_folder, pdb_id + '.mtz')
  cmd_list = []
  cmd_list.append('phenix.cif_as_mtz')
  cmd_list.append(cif)
  cmd_list.append('--output-file-name={}'.format(mtz_fn))
  cmd_list.append("--merge")
  cmd_list.append("--remove-systematic-absences")
  cmd_list.append("--map-to-asu")
  cmd_list.append("--ignore-bad-sigmas")
  cmd_list.append("--extend-flags")
  cmd = ' '.join(cmd_list)
  r = easy_run.go(cmd)
  # NOTE !!! in windows r does not returns the errors as expected
  tmp = [x for x in r.stdout_lines if '--' in x]
  tmp2 = ''.join(tmp)
  run_cmd_again = False
  if '--incompatible_flags_to_work_set' in tmp2:
    cmd_list.append('--incompatible_flags_to_work_set')
    run_cmd_again = True
  if '--symmetry' in tmp2:
    cmd_list.append('--symmetry={}'.format(pdb))
    run_cmd_again = True
  if run_cmd_again:
    cmd = ' '.join(cmd_list)
    easy_run.go(cmd)
  try:
    # Get miller arrays from mtz file
    mtz_object = iotbx.mtz.object(file_name=mtz_fn)
    miller_arrays = mtz_object.as_miller_arrays()
  except:
    miller_arrays = None
  # cleanup
  os.remove(cif)
  if not miller_arrays:
    return None
  return miller_arrays

def make_r_free_boolean(r_free_flags):
  '''
  Convert r_free_flag from any of the conventions possible
  to a boolean

  Posiblle convention for free and working set flags:
  CCP4     assigns the flag r_free_flags to be 0 for the free set and 1,
           ...n-1 for the working set.
  XPLOR    assigns the flag TEST to be 1 for the free set and 0 for the
           working set.
  CNS      assigns the flag TEST to be 1 for the free set and 0,2,...n-1 for
           the working set.
  SHELX    assigns a flag with -1 for the free set and 1 for the working set.
  TNT      assigns a flag with 0 to indicate the free set.
  '''
  flag_value = iotbx.reflection_file_utils.guess_r_free_flag_value(
    miller_array =r_free_flags)
  if flag_value is None:
    return None
  else:
    return r_free_flags.customized_copy(data=r_free_flags.data()==flag_value)

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
    elif (not (v is None)) and (v.__str__() != '0') and (v != []) and (v != {}):
      out_lst.append(template.format(k,v.__str__()))
  return out_lst
