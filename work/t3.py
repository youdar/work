from __future__ import division
from libtbx.utils import Sorry
from iotbx import pdb
from iotbx import cif
from cctbx import miller
from iotbx import crystal_symmetry_from_any
from iotbx import mtz
from iotbx.option_parser import option_parser
import os, sys
from libtbx import easy_run

#import cProfile

'''
 pdb_name = "1imh.pdb"
>>> pdb_inp = pdb.input(file_name=pdb_name)
>>> structure = pdb_inp.xray_structure_simple()
>>> f_miller = structure.structure_factors(d_min=2.85).f_calc()
'''


def run(sf_file_name,pdb_file_name):
  # check if files exist
  if not isfile(sf_file_name): raise Sorry('{} is not a file'.format(sf_file_name))
  if not isfile(pdb_file_name): raise Sorry('{} is not a file'.format(pdb_file_name))
  # start processing file
  # convert cif to mtz and process mtz file
  file_name_mtz = sf_file_name + '.mtz'
  # creates file in current folder
  easy_run.call("phenix.cif_as_mtz %s"%sf_file_name)
  mtz_object = mtz.object(file_name=file_name_mtz)
  # Process the mtz_object
  miller_arrays_dict = mtz_object.as_miller_arrays_dict()
  miller_arrays = mtz_object.as_miller_arrays()
  for x in miller_arrays_dict:
    print x
    print '+'*60
    #print miller_arrays_dict[x].show_array()
    print '='*60

  # read pdb info
  pdb_inp = pdb.input(file_name=pdb_file_name)          	# read the pdb file data
  structure = pdb_inp.xray_structure_simple()
  xray_structure_simple = pdb_inp.xray_structures_simple()	# a list of structure factors
  f_miller = xray_structure_simple[0].structure_factors(d_min=2.85).f_calc()

  # delete file
  os.remove(file_name_mtz)

  #


  print 'wait here'

def isfile(s):
  ''' (str)->boolean
  Check if s is a file in the current directory

  >>> isfile(s)
  '''
  return os.path.isfile(s)

def set_working_dir():
  #os.chdir('/net/cci-filer2/raid1/home/youval/Work/work')
  #os.chdir('/net/cci-filer2/raid1/home/youval/Work/work/from pavel/')
  os.chdir('c:\\Phenix\\Dev\Work\\work\\from pavel')
  print 'Current working directory:'
  print os.getcwd()
  print '='*50

if __name__=='__main__':
  set_working_dir()
  #sf_file_name = '/net/cci/pdb_mirror/structure_factors/00/r100dsf.ent.gz'
  #pdb_file_name = '/net/chevy/raid1/pdb_mirror/pdb/00/pdb100d.ent.gz'
  #run()
  #run(sf_file_name,pdb_file_name)
  #run('1akg-sf.cif.gz','pdb100d.ent.gz')
  #run('r1pglsf.ent.gz','pdb1pgl.ent.gz')
  run('r4iw4sf.ent.gz','pdb4iw4.ent.gz')
  #mtz_object = mtz.object(file_name=file_name)
  #return mtz_object.crystals()[0].crystal_symmetry()
  #print 'done'