from __future__ import division
from mmtbx.tls.tools import tls_from_pdb_inp
from iotbx import pdb
import os

__author__ = 'Youval'

"""
This file demonstrate how to scan all the pdb when working on
LBL machine. The scan is on the LBL pdb mirror
"""

def run():
  # set environment
  pdb_dir = os.environ["PDB_MIRROR_PDB"]
  pdb_files = open(os.path.join(pdb_dir, "INDEX"), "r").read().splitlines()
  # get all files including path
  pdb_files = [os.path.join(pdb_dir,f) for f in pdb_files]

  # To get information about the pdb...
  # example for one file
  fn = pdb_files[0]
  exp_type = get_experment_type
  # Get the 4 letter pdb code from the file name and path
  i = fn.find('pdb')
  pdb_code = fn[i+3:i+7]
  print 'Experiment type for {} is: {}'.format(pdb_code,exp_type)

def get_remark_iii_records(file_name,type_num):
  """ Collect remarks related to experiment type from pdb  """
  pdb_inp = pdb.hierarchy.input(file_name=file_name)
  pdb_inp_tls = tls_from_pdb_inp(
    remark_3_records = pdb_inp.input.extract_remark_iii_records(type_num),
    pdb_hierarchy = pdb_inp.hierarchy)
  return pdb_inp_tls.remark_3_records

def get_experment_type(file_name):
  '''(str) -> str
  Look for EXPERIMENT TYPE in PDB REMARK
  REMARK 200 		: X-RAY DIFFRACTION
  REMARK 210,215,217 	: NMR
  REMARK 230		: NEUTRON DIFFRACTION
  REMARK 245		: ELECTRON MICROSCOPE
  REMARK 250		: Other
  REMARK 265		: SMALL ANGLE X-RAY SCATTERING

  returns a string combined from all experiment types

  Note that the experiment identification is not done by actually reading the
  records, but by the record clasification at
  http://www.wwpdb.org/documentation/format33/remarks1.html
  '''
  remark_dict = dict(
    [(200,'X-RAY DIFFRACTION'),(210,'NMR'),(215,'NMR'),(217,'NMR'),
     (230,'NEUTRON DIFFRACTION'), (245,'ELECTRON MICROSCOPE'),
     (250,'Other'), (265,'SMALL ANGLE X-RAY SCATTERING')])
  experment_type = []
  for type_num in remark_dict.iterkeys():
    rem_records = get_remark_iii_records(file_name,type_num)
    if rem_records != []:
      experment_type.append(remark_dict[type_num])
  return experment_type

if __name__ == '__main__':
  run()
