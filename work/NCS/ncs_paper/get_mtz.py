from __future__ import division
import collect_ncs_files
import shutil
import sys
import os

__author__ = 'youva_000'

def run(pdb_id):
  """
  Make mtz file and updates ncs information records

  ASU need to exist in
  '/net/cci-filer2/raid1/home/youval/work/work/NCS/ncs_paper
  /ncs_paper_data_files/asu'

  the file log_"pdb id" need to exist in
  '/net/cci-filer2/raid1/home/youval/work/work/NCS/ncs_paper
  /ncs_queue_results'

  If there are no issues with the xray mtz file will be writen in
  '/net/cci-filer2/raid1/home/youval/work/work/NCS/ncs_paper
  /ncs_paper_data_files/mtz'

  if there are issues with the xray data, move pdb info records to
  '/net/cci-filer2/raid1/home/youval/work/work/NCS/ncs_paper
  /pdb_with_ncs_not_used'
  """
  assert len(pdb_id) == 1
  pdb_id = pdb_id[0]
  collect = collect_ncs_files.ncs_paper_data_collection()
  pdb = os.path.join(collect.asu_dir, pdb_id + '.pdb')
  pdb_info_fn = os.path.join(collect.pdb_records_dir,'log_' + pdb_id)
  pdb_info = collect.read_from_file(pdb_info_fn)
  pdb_info = fix_pdb_info(pdb_info)
  if pdb_info:
    pdb_info = collect.make_mtz_file(pdb_info)
    if pdb_info:
      collect.write_to_file('log_' + pdb_id,pdb_info)
    else:
      shutil.move(pdb_info_fn,collect.pdb_not_used_dir)

def fix_pdb_info(pdb_info):
  """ if using pdb_info record with old names, update them """

  pdb_info = collect_ncs_files.File_records()
  if hasattr(pdb_info, 'n_atoms_in_asu'):
    return pdb_info
  else:
    pdb_info_new = collect_ncs_files.File_records()
    for key in pdb_info.__dict__.iterkeys():
      if hasattr(pdb_info_new,key):
        pdb_info_new.__dict__[key] = pdb_info.__dict__[key]
      else:
        if key == 'solvent':
          new_key = 'solvent_fraction'
          pdb_info_new.__dict__[new_key] = pdb_info.__dict__[key]
        elif key == 'data_to_param_ration':
          new_key = 'data_to_param_ratio'
          pdb_info_new.__dict__[new_key] = pdb_info.__dict__[key]
        print 'Need to change the key: ',key
    return pdb_info_new


if __name__=='__main__':
  run(sys.argv[1:])
