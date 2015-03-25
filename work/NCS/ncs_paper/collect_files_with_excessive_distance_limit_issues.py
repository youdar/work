from __future__ import division
import collect_ncs_files
from libtbx.command_line import easy_qsub
from glob import glob
import shutil
import sys
import csv
import os

"""
Collect files that had issues with excessive_distance_limit
"""

def run():
  c = collect_ncs_files.ncs_paper_data_collection()
  fn = c.ncs_dir
  fn = os.path.join(fn,'ncs_paper_data.csv')
  with open(fn, 'rb') as csvfile:
    spamreader = csv.reader(csvfile, delimiter=' ', quotechar='|')
    for row in spamreader:
      print ', '.join(row)


  # refinement_dir_list = [
  #   'refine_no_ncs_dir','refine_cartesian_ncs','refine_torsion_ncs',
  #   'refine_ncs_con_no_oper','refine_ncs_con_all']
  # options = [
  #   '-no_ncs','-cartesian_ncs_restraints','-torsion_ncs_restraints',
  #   '-ncs_constraints_no_operators','-ncs_constraints_all']
  #
  # refinement_types = collect_ncs_files.get_refine_test_names()
  # #
  # pdb_ids = []
  # for rdn,rt,op in zip(refinement_dir_list,refinement_types,options):
  #   out_folder = c.__dict__[rdn]
  #   pdb_id_dirs = glob(os.path.join(out_folder,'*'))
  #   print len(pdb_id_dirs)
  #   for pdb_dir in pdb_id_dirs:
  #     pdb_id = pdb_dir[-4:]
  #     pdb_info_fn = os.path.join(c.data_dir,'log_' + pdb_id)
  #     pdb_info = c.read_from_file(pdb_info_fn)
  #     if pdb_info.refinement_records.has_key(rt):
  #       if pdb_info.refinement_records[rt].r_free_final == 0:
  #         # del: the line below will prevent redoing when -no_ncs has no value
  #         if pdb_info.refinement_records.has_key('no ncs'):
  #           if pdb_info.refinement_records['no ncs'].r_free_final > 0:
  #             # no result - rerun refinement
  #             pdb_ids.append(pdb_id)
  # print pdb_ids

if __name__ == '__main__':
  run()
