from __future__ import division
from compare_clashscores_CCTBX_vs_PROBE import get_pdb_id
import cPickle as pickle
import os

def run():
  """
  After running submit_clashscore_all_pdb_to_queue.py and
  collecting_and_looking_at_data.py some files might not have run.
  Collect the list of those files and create a text file containing that list

  The test file containing the PDB IDs of files that were not run is:
  "files_to_run.txt"

  PDB IDs of files that were processed are found in the pickled file "test_data"

  Note: This must be done on LBL machine, since we are looking at the index of
  PDB files on the LBL PDB mirror, "PDB_MIRROR_PDB"
  """
  path = '/net/cci-filer2/raid1/home/youval/Work/work/Clashes'
  fn = 'test_data'
  file_name = os.path.join(path,fn)
  if os.path.isfile(file_name):
    # get PDB IDs list from PDB_MIRROR_PDB
    pdb_dir = os.environ["PDB_MIRROR_PDB"]
    files = open(os.path.join(pdb_dir, "INDEX"), "r").read().splitlines()
    pdb_ids_in_mirror = {get_pdb_id(x) for x in files}

    # get processed file list from test_data
    data = pickle.load(open(file_name,'r'))
    # the first record in the data is the PDB id
    data = {x[0] for x in data}
    not_processed_pdbs = pdb_ids_in_mirror - data
    if not_processed_pdbs:
      print "number of unprocessed pdbs: ",len(not_processed_pdbs)
      print "saving PDB id list to files_to_run.txt"
      not_processed_pdbs = '\n'.join(not_processed_pdbs)
      out_file = os.path.join(path,'files_to_run.txt')
      open(out_file,'w').write(not_processed_pdbs)
    else:
      print "All files in LBL mirror were processed"
  else:
    print "test_data file does not exist"

if __name__ == '__main__':
  run()
