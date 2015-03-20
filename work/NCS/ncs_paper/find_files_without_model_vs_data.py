from __future__ import division
from datetime import datetime
import collect_ncs_files
from glob import glob
import cPickle as pickle
import unittest
import get_mtz
import shutil
import sys
import os


def run():
  c = collect_ncs_files.ncs_paper_data_collection()
  file_list = glob(os.path.join(c.model_vs_data_dir,'*.txt'))
  print "number of model vs data files:",len(file_list)
  empty_files = []
  for f in file_list:
    d = open(f,'r').read().splitlines()
    if len(d) < 5:
      empty_files.append(f[-8:-4])
  print empty_files
  print len(empty_files)


if __name__ == '__main__':
  run()
