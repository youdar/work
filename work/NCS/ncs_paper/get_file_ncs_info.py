from __future__ import division
import collect_ncs_files
import sys
import os

__author__ = 'Youval'


def run(pdb_id):
  """
  Collect information and saves ASU in the
  '/net/cci-filer2/raid1/home/youval/work/work/NCS/ncs_paper
  /ncs_paper_data_files/asu'

  Args:
    fn (pdb_id): pdb
  """
  assert len(pdb_id) == 1
  pdb_id = pdb_id[0]
  collect = collect_ncs_files.ncs_paper_data_collection()
  pdb_info = collect.get_pdb_file_info(pdb_id)
  if pdb_info:
    fn = os.path.join(collect.data_dir,'log_' + pdb_id)
    collect.write_to_file(fn,pdb_info)

if __name__=='__main__':
  run(sys.argv[1:])
