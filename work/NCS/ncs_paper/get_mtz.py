from __future__ import division
import collect_ncs_files
import shutil
import sys
import os

__author__ = 'Youval'

def run(pdb_id):
  """
  Make mtz file and updates ncs information records
  """
  assert len(pdb_id) == 1
  pdb_id = pdb_id[0]
  collect = collect_ncs_files.ncs_paper_data_collection()
  pdb_info_fn = os.path.join(collect.data_dir,'log_' + pdb_id)
  pdb_info = collect.read_from_file(pdb_info_fn)
  if pdb_info:
    pdb_info = collect.make_mtz_file(pdb_info)
    if pdb_info:
      collect.write_to_file(pdb_info_fn,pdb_info)
    else:
      shutil.move(pdb_info_fn,collect.pdb_not_used_dir)

if __name__=='__main__':
  run(sys.argv[1:])
