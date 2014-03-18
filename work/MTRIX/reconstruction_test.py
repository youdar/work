from __future__ import division
from  misc_scripts.r_factor_calc import *
import os


#os.chdir('/net/cci-filer2/raid1/home/youval/Work/work/')
os.chdir('/net/cci-filer2/raid1/home/youval/Work/work/MTRIX/')
print os.getcwd()
print r_factor_calc(['pdb4iw4.ent.gz','r4iw4sf.ent.gz'])
