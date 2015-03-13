from __future__ import division
from libtbx import easy_run
import collect_ncs_files
import sys
import os

__author__ = 'Youval'

def run(pdb_id):
  """  Collect info from phenix.ready_set  """
  assert len(pdb_id) == 1
  pdb_id = pdb_id[0]
  c = collect_ncs_files.ncs_paper_data_collection()
  pdb = os.path.join(c.asu_dir,pdb_id + '.pdb')
  cmd = "phenix.ready_set {} --silent".format(pdb)
  eff = os.path.join(c.cif_dir,pdb_id + '.eff')
  pdb_updated = eff.replace('.eff','.updated.pdb')
  current_dir = os.getcwd()
  os.chdir(c.cif_dir)
  r = easy_run.fully_buffered(cmd)
  if os.path.isfile(eff): os.remove(eff)
  if os.path.isfile(pdb_updated): os.remove(pdb_updated)
  os.chdir(current_dir)

if __name__=='__main__':
  run(sys.argv[1:])
