from __future__ import division
from libtbx import easy_run
import collect_ncs_files
import re
import shutil
import sys
import os

__author__ = 'Youval'

def run(pdb_id):
  """
  Collect info from mmtbx.model_vs_data
  """
  assert len(pdb_id) == 1
  pdb_id = pdb_id[0]
  c = collect_ncs_files.ncs_paper_data_collection()
  pdb_info_fn = os.path.join(c.data_dir,'log_' + pdb_id)
  pdb_info = c.read_from_file(pdb_info_fn)
  pdb = os.path.join(c.asu_dir,pdb_id + '.pdb')
  mtz = os.path.join(c.mtz_dir,pdb_id + '.mtz')
  out_file_name = os.path.join(c.model_vs_data_dir,pdb_id + '.txt')
  cmd = 'mmtbx.model_vs_data {} {} > {}'.format(pdb,mtz,out_file_name)
  r = easy_run.go(cmd)
  if 'Multiple equally' in r.stdout_lines[0]:
    cmd = 'mmtbx.model_vs_data {} {} f_obs_label=" F(+)" > {}'
    cmd.format(pdb,mtz,out_file_name)
    r = easy_run.go(cmd)
  #
  d = open(out_file_name,'r').read().splitlines()
  for l in d:
    match_work = re.search(r'(r_work.*re-computed.*:)(.*)',l)
    match_free = re.search(r'(r_free.*re-computed.*:)(.*)',l)
    if match_work:
      pdb_info.r_work_model_vs_data = float(match_work.group(2))
    if match_free:
      pdb_info.r_free_model_vs_data = float(match_free.group(2))
  # round some values
  pdb_info.data_completeness = round(pdb_info.data_completeness,3)
  pdb_info.data_to_param_ratio = round(pdb_info.data_to_param_ratio,3)
  pdb_info.solvent_fraction = round(pdb_info.solvent_fraction,3)
  # save the data collected
  c.write_to_file(pdb_info_fn,pdb_info)

if __name__=='__main__':
  run(sys.argv[1:])
