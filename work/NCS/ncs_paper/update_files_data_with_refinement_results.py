from __future__ import division
import collect_ncs_files
import sys

__author__ = 'Youval'

def run(args):
  """
  update a single or all PDB structure data, by collecting refinements
  results
  """
  pdb_id = None
  if len(args) == 1:
    pdb_id = args[0]
  c = collect_ncs_files.ncs_paper_data_collection()
  c.collect_refinement_results(pdb_id)
  print 'Done...'

if __name__ == '__main__':
  run(sys.argv[1:])
