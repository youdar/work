from __future__ import division
from libtbx.command_line import easy_qsub
from iotbx import pdb
import cProfile
import pstats
import iotbx.ncs
import time
import sys
import os

class null_out(object):
  """Pseudo-filehandle for suppressing printed output."""
  def isatty(self): return False
  def close(self): pass
  def flush(self): pass
  def write(self, str): pass
  def writelines(self, sequence): pass

def check_for_ncs(pdb_code=''):
  """
  when pdb_code='':
    Scan all pdb for structure containing NCS relations
  Otherwise process a single pdb file
  """
  assert isinstance(pdb_code,str)
  if (not pdb_code) or (pdb_code.lower() == 'all'):
    # Run on all PDB - Only on LBL machine
    osType = sys.platform
    msg = 'Please run this only on LBL computer'
    assert not osType.startswith('win'),msg
    # set environment
    pdb_dir = os.environ["PDB_MIRROR_PDB"]
    pdb_files = open(os.path.join(pdb_dir, "INDEX"), "r").read().splitlines()
    # send pdb_files to queue
    pdb_files = [os.path.join(pdb_dir,f) for f in pdb_files]
    submit_files_to_queue(pdb_files)
  else:
    # Process a single file
    pdb_dir = ''
    fn = get_pdb_file(file_name=pdb_code,print_out=False)
    pdb_code, n_mtrix_rec, d = check_a_single_file(fn)
    # d:
    # t_min, t_max, t_simple,
    # t_min_gruops, t_min_gruops, t_simple_gruops,
    # t_min_time, t_max_time, t_simple_time
    s1 =  '{0:<6}   {1:<6}   {2:<6}  {3:<6}'
    s2 =  '{0:<6.3f}  {1:<6.3f}  {2:<6.3f}'
    s3 =  '{0:<6}   {1:<6}   {2:<6} '
    print 'Results for %s are:'%pdb_code
    print s1.format('min','max','simple','in pdb')
    print s1.format(d[0], d[1], d[2], n_mtrix_rec)
    print 'Number of NCS groups'
    print s3.format(d[3],d[4],d[5])
    print 'In time'
    print s2.format(d[6],d[7],d[8])
    # clean files
    if not ('pdb_mirror' in fn):
      if os.path.isfile(fn): os.remove(fn)

def get_pdb_file(file_name, print_out=True):
  """ (file_name) -> file_path
  This function will check if a pdb file_name exist.
  If it is not, it will either fetch it from the RCSB website
  or find it on LBLs pdb mirror folder

  Args:
    file_name (str): a pdb file name
  Returns:
    file_name (str): the location, path, of the file
  """
  if not os.path.isfile(file_name):
    # get a clean pdb file name
    if print_out:
      s = 'No such file in working directory. ' \
          'Trying to fetch {} file from RCSB web site'
      print s.format(file_name)
    file_name = get_4_letters_pdb_id(file_name)
    osType = sys.platform
    if osType.startswith('win'):
      from iotbx.pdb import fetch
      # fetch pdb file from intenet
      file_name = fetch.get_pdb (file_name,'pdb',mirror='rcsb',log=null_out())
    else:
      # find the file in LBL pdb mirror folder
      pdb_dir = os.environ["PDB_MIRROR_PDB"]
      pdb_files = open(os.path.join(pdb_dir, "INDEX"), "r").read().splitlines()
      for i,p in enumerate(pdb_files):
        if file_name in p:
          break
      file_name = os.path.join(pdb_dir,pdb_files[i])
      # # Second method
      # f = 'pdb{}.ent.gz'.format(file_name)
      # file_name = []
      # for root, _, files in os.walk(pdb_dir):
      #   if f in files:
      #     file_name = os.path.join(root, f)
      #     break
  elif print_out:
    print 'Using the file {} found in the working directory'.format(file_name)
  return file_name

def submit_files_to_queue(pdb_files):
  """
  process all files in pdb_files using the queuing system
  Output will be files with names "log_xxxx" in the current directory
  """
  # Set command path
  phenix_source = "/net/chevy/raid1/youval/Work_chevy/phenix_build/setpaths.csh"
  path = '/net/cci/youval/Work/work/NCS/pdb_surveys/'
  com_path = path + 'compare_find_ncs_to_mtrix_records.py'
  where_to_run_dir = os.getcwd()
  commands = []
  for fn in pdb_files:
    pdb_code = get_4_letters_pdb_id(fn)
    s = 'python {0} {1} >& log_{2}'.format(com_path,fn,pdb_code)
    commands.append(s)
  # Send to queue
  easy_qsub.run(
    phenix_source = phenix_source,
    where         = where_to_run_dir,
    commands      = commands,
    qsub_cmd      = 'qsub -q all.q@morse',
    size_of_chunks= 300)

def check_a_single_file(file_name_and_path):
  """  Process a single file

  Args:
    file_name_and_path (str): pdb file to process, including path

  Returns:
    data: t_min, t_max, t_simple, t_min_gruops, t_min_gruops, t_simple_gruops,
      t_min_time, t_max_time, t_simple_time
      t_min (int): min number of transforms found in pdb file
      t_max (int): max number of transforms found in pdb file
      t_simple (int): max number of transforms found by simple_ncs_from_pdb
      t_..._gruops (int): number of ncs groups by each method)
      t_..._time (float): time in sec, each method run

    n_mtrix_rec (int): the number of not-present MTRIX records in the PDB file
  """
  pdb_code = get_4_letters_pdb_id(file_name_and_path)
  n_mtrix_rec = number_of_mtrix_rec(file_name_and_path)
  # Check for NCS operators
  data = test(file_name_and_path)
  return pdb_code, n_mtrix_rec, data

def number_of_mtrix_rec(file_name_and_path,eps=0.01):
  """  Get the number of MTRIX records in PDB file
  Args:
    file_name_and_path (str)
    eps (float): Allowed error when testing rotation matrices

  Returns:
    file_name_and_path (str): file name and path of the complete ASU
    n_mtrix_rec (int) the number of not-present MTRIX records in the PDB file
"""
  pdb_hierarchy_inp = pdb.hierarchy.input(file_name=file_name_and_path)
  tr_info = pdb_hierarchy_inp.input.process_mtrix_records(eps=eps)
  n_mtrix_rec = 0
  for r,t in zip(tr_info.r, tr_info.t):
    if not (r.is_r3_identity_matrix() and t.is_col_zero()):
      n_mtrix_rec += 1
  return n_mtrix_rec

def get_4_letters_pdb_id(file_name):
  """(str)  -> str
  clean a pdb file name, remove path and file extensions

  :param file_name (str): pdb file name that may look like pdb1a37.pdb
  :return pdb_id (str): the 4 letter pdb id

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

def test(file_name_and_path):
  """
  Check if NCS relations exist and return the minimal and maximal number of
  transforms

  Args:
    file_name_and_path (str): file name, including path

  Returns:
    min_ncs_opr, max_ncs_opr (int, int): min and max number of transforms
  """
  t_min = t_max = t_simple = -1
  t_min_groups = t_max_groups = t_simple_groups = -1

  time_min = time.time()
  try:
    trans_obj = iotbx.ncs.input(
        file_name=file_name_and_path,
        use_minimal_master_ncs=False,
        write_messages=True,
        process_similar_chains=True,
        min_percent=0.85,
        max_rmsd=99999)
    t_min = len(trans_obj.transform_to_ncs)
    t_min_groups = trans_obj.number_of_ncs_groups
    trans_obj.get_ncs_info_as_spec(write=True)
    print '='*30
  except:
    pass
  t_min_time = (time.time()-time_min)

  time_max = time.time()
  # try:
  trans_obj = iotbx.ncs.input(
      file_name=file_name_and_path,
      use_minimal_master_ncs=True,
      write_messages=False,
      process_similar_chains=True,
      min_percent=0.85,
      max_rmsd=99999)
  t_max = len(trans_obj.transform_to_ncs)
  t_max_groups = trans_obj.number_of_ncs_groups
  trans_obj.get_ncs_info_as_spec(write=True)
  print '='*30
  # except:
  #   pass
  t_max_time = (time.time()-time_max)

  # Fixme: change use_simple_ncs_from_pdb back to False
  time_simple = time.time()
  try:
    trans_obj = iotbx.ncs.input(
        file_name=file_name_and_path,
        write_messages=False,
        process_similar_chains=False,
        min_percent=0.85,
        max_rmsd=99999)
    t_simple = len(trans_obj.transform_to_ncs)
    t_simple_groups = trans_obj.number_of_ncs_groups
    trans_obj.get_ncs_info_as_spec(write=True)
    print '='*30
  except:
    pass
  t_simple_time = (time.time()-time_simple)

  return t_min, t_max, t_simple, \
         t_min_groups, t_max_groups, t_simple_groups,\
         t_min_time, t_max_time, t_simple_time

def run(args):
  """
  Run an all pdb test
  1) Work on a LBL machine
  2) use check_for_ncs(pdb_code='') if the function call below

  You can run a test on a single pdb file on any machine
  """
  if len(args) == 0:
    if sys.platform.startswith('win'):
      dr = r'C:\Phenix\Dev\Work\work\NCS\pdb_surveys\pdb_survey_logs'
    else:
      dr = '/net/cci-filer2/raid1/home/youval/Work/work/NCS/pdb_surveys/pdb_survey_logs'
    os.chdir(dr)
    pc=''
    # pc='2a2u'
    # pc='4bed'
    # pc='1k5j'
    # pc='4erm'
    # pc='1a5d'
    # pc='4b31'
    # pc='1k5d'
    # pc='200d'
    # pc='2g32'
    # pc='4hs6'
    # pc='1jj2'
    # pc='1ruy'
    # pc='4d0i'
    # pc='1z7q'
    # pc='3abh'
    # pc='4bnz'
    # pc='4bq9'
    pc='4ang'
    if pc:
      check_for_ncs(pdb_code=pc)
      # s = 'check_for_ncs(pdb_code="%s")' % pc
      # cProfile.run(s,filename='cProfile_log')
      # p = pstats.Stats('cProfile_log')
      # p.sort_stats('time').print_stats(15)
    else:
      # for all pdb test use check_for_ncs(pdb_code=pc)
      check_for_ncs(pdb_code='')
  elif len(args) == 1:
    fn = args[0]
    pdb_code, n_mtrix_rec, out_data = check_a_single_file(fn)
    # out_data
    # t_min, t_max, t_simple,
    # t_min_gruops, t_min_gruops, t_simple_gruops,
    # t_min_time, t_max_time, t_simple_time
    s = '{0},{1},{2},{3},{4},{5},{6},{7},{8:0.3f},{9:0.3f},{10:0.3f}'
    print s.format(pdb_code, n_mtrix_rec, *out_data)

if __name__ == '__main__':
  run(sys.argv[1:])

