import os
from iotbx import pdb

def run():
  os.chdir(r'c:\Phenix\Dev\Work\work\junk')
  pdb_inp = pdb.input(file_name='4dzi.pdb')
  pdb_inp.process_BIOMT_records()

if __name__=='__main__':
  run()