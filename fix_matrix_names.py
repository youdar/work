from __future__ import division
from  misc_scripts.r_factor_calc import *
from  iotbx.pdb.multimer_reconstruction import multimer
from iotbx import pdb
import cPickle as pickle
import os


os.chdir('/net/cci-filer2/raid1/home/youval/Work/work')
MTRIX_with_Straucture_Factor = pickle.load(open('MTRIX_with_Straucture_Factor_file_list','r'))

s = set()
for x in MTRIX_with_Straucture_Factor:
    if '/' not in x:
        s.add(x)

pickle.dump(list(s),open('MTRIX_with_Straucture_Factor_file_list','w'))