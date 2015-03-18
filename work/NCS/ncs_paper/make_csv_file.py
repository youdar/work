from __future__ import division
import collect_ncs_files
import sys
import os

__author__ = 'Youval'

def run(args):
  """
  Collect data on all tested PDB structures and save it in
  make_csv_file.py

  Args
    out_path (str): alternative path for output file
  """
  assert len(args) < 2
  file_name = 'ncs_paper_data.csv'
  c = collect_ncs_files.ncs_paper_data_collection()
  if args:
    out_path = args[0]
    fn = os.path.join(out_path,file_name)
  else:
    out_path = ''
    fn = os.path.join(c.ncs_dir,file_name)
  c.make_csv_file(file_name=file_name,out_path=out_path)
  print 'Save csv file to\n{}\n'.format(fn)

if __name__ == '__main__':
  run(sys.argv[1:])
