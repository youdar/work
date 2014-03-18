import os
import iotbx.pdb

def full_pdb_file_paths():
  pdb_dir = os.environ["PDB_MIRROR_PDB"]
  pdb_files = open(os.path.join(pdb_dir, "INDEX"), "r").readlines()
  result = []
  for pdb_file in pdb_files:
    pdb_code = pdb_file[-12:-8]
    pdb_file = os.path.join(pdb_dir, pdb_file.strip())
    result.append(pdb_file)
  return result

def run():
  for i_file, file_name in enumerate(full_pdb_file_paths()):
    print "processing:", i_file, file_name
    try: pdb_inp = iotbx.pdb.input(file_name = file_name)
    except Exception, e: 
      print "FAILED: %s"%str(e)
    biomt_obj = pdb_inp.process_BIOMT_records()
    mtrix_obj = pdb_inp.process_mtrix_records()

if (__name__ == "__main__"):
  run()
