from __future__ import division
import collect_ncs_files
import sys
import os

__author__ = 'Youval'


def run(args):
  """
  Print information collected on PDB files

  Usage example:
  >>> python see_pdb_info.py 1vcr 1pgf
  """
  if len(args) == 0:
    print run.__doc__
  else:
    collect = collect_ncs_files.ncs_paper_data_collection()
    for pdb_id in args:
      pdb_id = pdb_id.replace('.pdb','')
      pdb_info_fn = os.path.join(collect.data_dir,'log_' + pdb_id)
      if os.path.isfile(pdb_info_fn):
        pdb_info = collect.read_from_file(pdb_info_fn)
        print '                 ---    {}   ---'.format(pdb_id)
        print pdb_info
        print '.'*40
      else:
        print 'Could not find info for:',pdb_id

if __name__=='__main__':
  run(sys.argv[1:])
