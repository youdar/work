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
import re

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
    # model vs data
    self.r_free_model_vs_data = None
    self.r_work_model_vs_data = None

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
    self.clashscore_final = None
    self.c_beta_deviation = None
    self.c_beta_final = None
    self.map_cc = None
    self.rmsd = None
    self.rama_outliers = None
    self.rama_final = None
    self.rotamer_outliers = None
    self.rotamer_final = None

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
    self.cif_dir = os.path.join(self.ncs_dir,'cif')
    self.data_dir = os.path.join(self.ncs_dir,'data')
    self.refine_no_ncs_dir = os.path.join(self.ncs_dir,'refine_no_ncs')
    self.refine_cartesian_ncs = os.path.join(self.ncs_dir,'refine_cartesian_ncs')
    self.refine_torsion_ncs = os.path.join(self.ncs_dir,'refine_torsion_ncs')
    self.refine_ncs_con_no_oper = os.path.join(self.ncs_dir,'refine_ncs_no_oper')
    self.refine_ncs_con_all = os.path.join(self.ncs_dir,'refine_ncs_con_all')
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

    Args:
      file_name (str): such as log_"pdb id"
      file_record (obj): file record object
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
    files_list = glob(os.path.join(self.data_dir,'log_*'))
    for fn in files_list:
      r = pickle.load(open(fn,'r'))
      self.pdbs_dict[r.pdb_id] = r
    self.files_list = [x[-4:] for x in files_list]
    return self.pdbs_dict

  def collect_refinement_results(self):
    """ updates records with refinement results """
    paths = [
      self.refine_no_ncs_dir,self.refine_cartesian_ncs,
      self.refine_torsion_ncs,self.refine_ncs_con_no_oper,
      self.refine_ncs_con_all]
    refine_test_names = [
    'no ncs','cartesian ncs restraints','torsion ncs restraints',
    'ncs constraints no operators','ncs constraints all']
    records = self.collect_all_file_records()
    for test_name,data_path in zip(refine_test_names,paths):
      # get all folders in directory
      if os.path.isdir(data_path):
        pdb_id_dirs = glob(os.path.join(data_path,'*'))
        for pdb_dir in pdb_id_dirs:
          # update relevant record with new data
          pdb_id = pdb_dir[-4:]
          pdb_info = records[pdb_id]
          refine_results = collect_refine_data(pdb_dir)
          pdb_info.refinement_records[test_name] = refine_results
          fn = os.path.join(self.data_dir,'log_' + pdb_id)
          self.write_to_file(fn,pdb_info)

  def make_csv_file(self,file_name='',records=None,out_path=''):
    """
    creates a csv file from all data collected

    Args:
      file_name (str): output file name
      records (dict): dictionary containing all records
      out_path (str): output file path
      """
    # get pdb IDs to collect data on
    if not records:
      records = self.collect_all_file_records()
    if not records: return False
    if not file_name:
      file_name = 'ncs_paper_data.csv'
    if not out_path:
      out_path = self.ncs_dir
    headers, table_pos_map = table_headers()
    l = len(table_pos_map)
    h = [(v,k) for k,v in table_pos_map.iteritems()]
    h.sort()
    table = [[k for (v,k) in h]]
    file_rec = File_records()
    for pdb_id in records:
      pdb_info = records[pdb_id]
      new_row = ['',] * l
      for key in file_rec.__dict__.iterkeys():
        # fixme: remove the following test since all keys should be present
        if  pdb_info.__dict__.has_key(key):
          v = pdb_info.__dict__[key]
        else:
          v = None
        if not (v is None):
          if key == 'refinement_records':
            # unpack the dictionary containing different refinement tests
            for ref_type in v.iterkeys():
              ref_rec = v[ref_type]
              if not (ref_rec is None):
                # iterate over refinement test results
                for ref_k in ref_rec.__dict__.iterkeys():
                    head = headers.refinement_records[ref_type]
                    h = head.__dict__[ref_k]
                    if not (h is None):
                      i = table_pos_map[h]
                      d = ref_rec.__dict__[ref_k]
                      if d is None: d = ''
                      new_row[i] = str(d)
          elif key != 'issues':
            h = headers.__dict__[key]
            i = table_pos_map[h]
            new_row[i] = str(v)
      table.append(new_row)
    table = [','.join(x) for x in table]
    table = '\n'.join(table)
    fn = os.path.join(out_path,file_name)
    open(fn,'w').write(table)
    return True

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

def get_dict_as_list(d,template,add_title=False):
  """
  recursively expands dictionary for printing

  Args:
    d (dict):  a dictionary
    template (str): a template used to format the key and value of dict
    add_title (bool): add title to dictionary printout
  Returns:
    out_lst (list): a list of string containing the formatted key-value pairs
  """
  out_lst = []
  keys = d.keys()
  keys.sort()
  for k in keys:
    if add_title:
      out_lst.extend([k,'-'*len(k)])
    v = d[k]
    if type(v) is dict:
      x = get_dict_as_list(v,template,add_title=True)
      out_lst.extend(x)
    elif 'refinement_records' in type(v).__name__:
      out_lst.append(v.__repr__())
    elif (not (v is None)) and (v.__str__() != '0') and (v != []) and (v != {}):
      out_lst.append(template.format(k,v.__str__()))
    if add_title:
      out_lst.append('-'*40)
  return out_lst

def collect_refine_data(pdb_dir):
  """ collecting data from refinement log """
  refine_results = Refinement_results()
  files_list = glob(os.path.join(pdb_dir,'*.log'))
  if len(files_list) > 1:
    msg = "There are several refinement log files in: \n{}\n"
    msg += "please remove the .log extension from the files you do not collect"
    print msg.format(pdb_dir)
  elif files_list:
    data = open(files_list[0],'r').read().splitlines()
    i = 0
    for l in data:
      # collect data from file
      start_r_val = re.search(r'(Start R-work =)(.*)(\,.*R-free =)(.*)',l)
      final_r_val = re.search(r'(Final R-work =)(.*)(\,.*R-free =)(.*)',l)
      clashscore = re.search(r'(all-atom clashscore.*:)(.*)',l)
      rotamer_outliers = re.search(r'(rotamer outliers.*:)(.*)(\%)',l)
      c_beta_deviation = re.search(r'(cbeta deviations.*:)(.*)',l)
      cpu_time = re.search(r'(Total CPU time:)(.*)(minutes)',l)
      next_line_are_results = re.search(r'Accepted refinement result:',l)
      molprobity_statistics = re.search(r'Molprobity statistics',l)
      i += 1
      # update record
      if start_r_val:
        refine_results.r_work_init = float(start_r_val.group(2))
        refine_results.r_free_init = float(start_r_val.group(4))
      if final_r_val:
        refine_results.r_work_final = float(final_r_val.group(2))
        refine_results.r_free_final = float(final_r_val.group(4))
      if cpu_time:
        # convert cpu time to seconds
        refine_results.refinement_time = round(60*float(cpu_time.group(2)),1)
      if clashscore:
        refine_results.clashscore = float(clashscore.group(2))
      if rotamer_outliers:
        refine_results.rotamer_outliers = float(rotamer_outliers.group(2))
      if c_beta_deviation:
        refine_results.c_beta_deviation = float(c_beta_deviation.group(2))
      if next_line_are_results:
        d = data[i].split()
        if len(d) == 12:
          refine_results.clashscore_final = float(d[5])
          refine_results.rama_final = float(d[6])
          refine_results.rotamer_final = float(d[7])
          refine_results.c_beta_final = float(d[8])
      if molprobity_statistics:
        d = data[i + 2].split()
        if d[0].lower() == 'outliers':
          refine_results.rama_outliers = float(d[2])
    return refine_results
  else:
    return None

def table_headers():
  """
  Returns:
    headers (dict): map data and experiment to header
    table_pos_map (dict): map header to a location in the table
  """
  headers = File_records()
  headers.pdb_id = 'pdb id'
  headers.n_ncs_copies = 'n copies'
  headers.year = 'year'
  headers.resolution =  'resolution'
  headers.data_completeness =  'completeness'
  headers.solvent_fraction =  'solvent fraction'
  headers.experiment_type =  'experiment'
  headers.only_master_in_pdb =  'master only'
  headers.n_atoms_in_asu =  'atoms in asu'
  headers.data_to_param_ratio_ncs =  'p-to-d ratio ncs'
  headers.data_to_param_ratio =  'p-to-d ratio asu'
  headers.r_free_header =  'r-free header'
  headers.r_work_header =  'r-work header'
  headers.r_free_model_vs_data = 'r-free model vs data'
  headers.r_work_model_vs_data = 'r-work model vs data'
  #
  test = Refinement_results()
  test.r_free_init        = 'r-free init : no ncs'
  test.r_work_init        = 'r-work init : no ncs'
  test.r_free_final       = 'r-free final : no ncs'
  test.r_work_final       = 'r-work final : no ncs'
  test.refinement_time    = 'refinement time : no ncs'
  test.clashscore         = 'all-atom clashscore : no ncs'
  test.clashscore_final   = 'final clashscore : no ncs'
  test.rotamer_outliers   = 'rotamer outliers : no ncs'
  test.rotamer_final      = 'rotamer final : no ncs'
  test.c_beta_deviation   = 'cbeta deviations : no ncs'
  test.c_beta_final       = 'cbeta final : no ncs'
  test.rama_outliers      = 'rama outliers : no ncs'
  test.rama_final         = 'rama final : no ncs'
  headers.refinement_records['no ncs'] = test
  #
  test = Refinement_results()
  test.r_free_init        = 'r-free init : cartesian ncs restraints'
  test.r_work_init        = 'r-work init : cartesian ncs restraints'
  test.r_free_final       = 'r-free final : cartesian ncs restraints'
  test.r_work_final       = 'r-work final : cartesian ncs restraints'
  test.refinement_time    = 'refinement time : cartesian ncs restraints'
  test.clashscore         = 'all-atom clashscore : cartesian ncs restraints'
  test.clashscore_final   = 'final clashscore : cartesian ncs restraints'
  test.rotamer_outliers   = 'rotamer outliers : cartesian ncs restraints'
  test.rotamer_final      = 'rotamer final : cartesian ncs restraints'
  test.c_beta_deviation   = 'cbeta deviations : cartesian ncs restraints'
  test.c_beta_final       = 'cbeta final : cartesian ncs restraints'
  test.rama_outliers      = 'rama outliers : cartesian ncs restraints'
  test.rama_final         = 'rama final : cartesian ncs restraints'
  headers.refinement_records['cartesian ncs restraints'] = test
  #
  test = Refinement_results()
  test.r_free_init        = 'r-free init : torsion ncs restraints'
  test.r_work_init        = 'r-work init : torsion ncs restraints'
  test.r_free_final       = 'r-free final : torsion ncs restraints'
  test.r_work_final       = 'r-work final : torsion ncs restraints'
  test.refinement_time    = 'refinement time : torsion ncs restraints'
  test.clashscore         = 'all-atom clashscore : torsion ncs restraints'
  test.clashscore_final   = 'final clashscore : torsion ncs restraints'
  test.rotamer_outliers   = 'rotamer outliers : torsion ncs restraints'
  test.rotamer_final      = 'rotamer final : torsion ncs restraints'
  test.c_beta_deviation   = 'cbeta deviations : torsion ncs restraints'
  test.c_beta_final       = 'cbeta final : torsion ncs restraints'
  test.rama_outliers      = 'rama outliers : torsion ncs restraints'
  test.rama_final         = 'rama final : torsion ncs restraints'
  headers.refinement_records['torsion ncs restraints'] = test
  #
  test = Refinement_results()
  test.r_free_init      = 'r-free init : ncs constraints no operators'
  test.r_work_init      = 'r-work init : ncs constraints no operators'
  test.r_free_final     = 'r-free final : ncs constraints no operators'
  test.r_work_final     = 'r-work final : ncs constraints no operators'
  test.refinement_time  = 'refinement time : ncs constraints no operators'
  test.clashscore       = 'all-atom clashscore : ncs constraints no operators'
  test.clashscore_final = 'final clashscore : ncs constraints no operators'
  test.rotamer_outliers = 'rotamer outliers : ncs constraints no operators'
  test.rotamer_final    = 'rotamer final : ncs constraints no operators'
  test.c_beta_deviation = 'cbeta deviations : ncs constraints no operators'
  test.c_beta_final     = 'cbeta final : ncs constraints no operators'
  test.rama_outliers    = 'rama outliers : ncs constraints no operators'
  test.rama_final       = 'rama final : ncs constraints no operators'
  headers.refinement_records['ncs constraints no operators'] = test
  #
  test = Refinement_results()
  test.r_free_init      = 'r-free init : ncs constraints all'
  test.r_work_init      = 'r-work init : ncs constraints all'
  test.r_free_final     = 'r-free final : ncs constraints all'
  test.r_work_final     = 'r-work final : ncs constraints all'
  test.refinement_time  = 'refinement time : ncs constraints all'
  test.clashscore       = 'all-atom clashscore : ncs constraints all'
  test.clashscore_final = 'final clashscore : ncs constraints all'
  test.rotamer_outliers = 'rotamer outliers : ncs constraints all'
  test.rotamer_final    = 'rotamer final : ncs constraints all'
  test.c_beta_deviation = 'cbeta deviations : ncs constraints all'
  test.c_beta_final     = 'cbeta final : ncs constraints all'
  test.rama_outliers    = 'rama outliers : ncs constraints all'
  test.rama_final       = 'rama final : ncs constraints all'
  headers.refinement_records['ncs constraints all'] = test
  #
  # map dictionary record to a location in the table
  headers_list = [
      'pdb id',
      'n copies',
      'year',
      'resolution',
      'completeness',
      'solvent fraction',
      'experiment',
      'master only',
      'atoms in asu',
      'p-to-d ratio ncs',
      'p-to-d ratio asu',
      'r-free header',
      'r-work header',
      'r-free model vs data',
      'r-work model vs data',
      'r-free init : no ncs',
      'r-work init : no ncs',
      'r-free final : no ncs',
      'r-work final : no ncs',
      'refinement time : no ncs',
      'all-atom clashscore : no ncs',
      'final clashscore : no ncs',
      'rotamer outliers : no ncs',
      'rotamer final : no ncs',
      'cbeta deviations : no ncs',
      'cbeta final : no ncs',
      'rama outliers : no ncs',
      'rama final : no ncs',
      'r-free init : cartesian ncs restraints',
      'r-work init : cartesian ncs restraints',
      'r-free final : cartesian ncs restraints',
      'r-work final : cartesian ncs restraints',
      'refinement time : cartesian ncs restraints',
      'all-atom clashscore : cartesian ncs restraints',
      'final clashscore : cartesian ncs restraints',
      'rotamer outliers : cartesian ncs restraints',
      'rotamer final : cartesian ncs restraints',
      'cbeta deviations : cartesian ncs restraints',
      'cbeta final : cartesian ncs restraints',
      'rama outliers : cartesian ncs restraints',
      'rama final : cartesian ncs restraints',
      'r-free init : torsion ncs restraints',
      'r-work init : torsion ncs restraints',
      'r-free final : torsion ncs restraints',
      'r-work final : torsion ncs restraints',
      'refinement time : torsion ncs restraints',
      'all-atom clashscore : torsion ncs restraints',
      'final clashscore : torsion ncs restraints',
      'rotamer outliers : torsion ncs restraints',
      'rotamer final : torsion ncs restraints',
      'cbeta deviations : torsion ncs restraints',
      'cbeta final : torsion ncs restraints',
      'rama outliers : torsion ncs restraints',
      'rama final : torsion ncs restraints',
      'r-free init : ncs constraints no operators',
      'r-work init : ncs constraints no operators',
      'r-free final : ncs constraints no operators',
      'r-work final : ncs constraints no operators',
      'refinement time : ncs constraints no operators',
      'all-atom clashscore : ncs constraints no operators',
      'final clashscore : ncs constraints no operators',
      'rotamer outliers : ncs constraints no operators',
      'rotamer final : ncs constraints no operators',
      'cbeta deviations : ncs constraints no operators',
      'cbeta final : ncs constraints no operators',
      'rama outliers : ncs constraints no operators',
      'rama final : ncs constraints no operators',
      'r-free init : ncs constraints all',
      'r-work init : ncs constraints all',
      'r-free final : ncs constraints all',
      'r-work final : ncs constraints all',
      'refinement time : ncs constraints all',
      'all-atom clashscore : ncs constraints all',
      'final clashscore : ncs constraints all',
      'rotamer outliers : ncs constraints all',
      'rotamer final : ncs constraints all',
      'cbeta deviations : ncs constraints all',
      'cbeta final : ncs constraints all',
      'rama outliers : ncs constraints all',
      'rama final : ncs constraints all',
    ]
  table_pos_map = {x:i for i,x in enumerate(headers_list)}
  return headers, table_pos_map


if __name__ == '__main__':
  """ update all PDB structure data, buy collecting refinements results """
  c = ncs_paper_data_collection()
  c.collect_refinement_results()
  print 'Done...'
