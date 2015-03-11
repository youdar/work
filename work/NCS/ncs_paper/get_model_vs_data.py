from __future__ import division
from libtbx import easy_run
import collect_ncs_files
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
  mtz = os.path.join(c.asu_dir,pdb_id + '.mtz')
  out_file_name = os.path.join(c.model_vs_data_dir,pdb_id + '.txt')
  cmd = 'mmtbx.model_vs_data {} {} > {}'.format(pdb,mtz,out_file_name)
  r = easy_run.go(cmd)
  tmp = [x for x in r.stdout_lines]
  print 'Done'

if __name__=='__main__':
  run(sys.argv[1:])
