import iotbx.cif
from cctbx import miller
from iotbx import crystal_symmetry_from_any
import sys,os

def run(args):
  cif_file, pdb_file = args
  print pdb_file
  print '*'*50
  cs = crystal_symmetry_from_any.extract_from(pdb_file)
  base_array_info = miller.array_info(
    crystal_symmetry_from_file=cs)
  all_miller_arrays = iotbx.cif.reader(file_path=cif_file).build_miller_arrays(
    base_array_info=base_array_info)
  print all_miller_arrays
  print 'wait here'

if (__name__ == "__main__"):
  os.chdir('/net/cci-filer2/raid1/home/youval/Work/work/from pavel/')
  run(sys.argv[1:])