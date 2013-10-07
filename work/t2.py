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


def run(sf_file_name,pdb_file_name):
  # check if files exist
  if not isfile(sf_file_name): raise Sorry('{} is not a file'.format(sf_file_name))
  if not isfile(pdb_file_name): raise Sorry('{} is not a file'.format(pdb_file_name))
  # start processing file
  cs = crystal_symmetry_from_any.extract_from(pdb_file_name)
  base_array_info = miller.array_info(crystal_symmetry_from_file=cs)
  all_miller_arrays = cif.reader(file_path=sf_file_name).build_miller_arrays(base_array_info=base_array_info)
  #
  for (data_name, miller_arrays) in all_miller_arrays.iteritems():
    print data_name
    for ma in miller_arrays.values():
      print get_label(ma),ma

  print 'wait here'

def isfile(s):
  ''' (str)->boolean
  Check if s is a file in the current directory

  >>> isfile(s)
  '''
  return os.path.isfile(s)

def get_label(miller_array):
  label = None
  for l in miller_array.info().labels:
    if ('_meas' in l) :
      if miller_array.is_xray_amplitude_array():
        label = "FOBS"
      elif miller_array.is_xray_intensity_array():
        label = "IOBS"
      break
    elif miller_array.anomalous_flag():
      if miller_array.is_xray_amplitude_array():
        label = "F"
      elif miller_array.is_xray_intensity_array():
        label = "I"
      break
    elif 'status' in l or '_free' in l:
      label = output_r_free_label
      break
    elif miller_array.is_hendrickson_lattman_array():
      label = "HL"
    elif (miller_array.is_complex_array()) :
      if (l.endswith("DELFWT")) :
        label = "DELFWT"
        break
      elif (l.endswith("FWT")) :
        label = "FWT"
        break
    elif (miller_array.is_real_array()) :
      if (l.endswith(".phase_calc")) :
        label = "PHIC"
        break
      elif ("pdbx_anom_difference" in l) :
        label = "DANO"
        break
  #if label is not None:
    #label_base = label
    #i = 1
    #while label in column_labels:
      #label = label_base + "-%i" %(i)
      #i += 1
  return label

def set_working_dir():
  #os.chdir('/net/cci-filer2/raid1/home/youval/Work/work')
  os.chdir('/net/cci-filer2/raid1/home/youval/Work/work/from pavel/')
  print 'Current working directory:'
  print os.getcwd()
  print '='*50

if __name__=='__main__':
  set_working_dir()
  #sf_file_name = '/net/cci/pdb_mirror/structure_factors/00/r100dsf.ent.gz'
  #pdb_file_name = '/net/chevy/raid1/pdb_mirror/pdb/00/pdb100d.ent.gz'

  #run(sf_file_name,pdb_file_name)
  #run('1akg-sf.cif.gz','pdb100d.ent.gz')
  #mtz_object = mtz.object(file_name=file_name)
  easy_run.call("phenix.cif_as_mtz %s"%'1akg-sf.cif.gz')
  mtz_object = mtz.object(file_name='1akg-sf.cif.gz.mtz')
  print 'done'
  #return mtz_object.crystals()[0].crystal_symmetry()