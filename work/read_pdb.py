import os, sys
from iotbx import pdb

def run(*args):
  print "run",args
  pdb_inp = pdb.input(args[0])
  print dir(pdb_inp)
  hierarchy = pdb_inp.construct_hierarchy()
  print dir(hierarchy)

if __name__=="__main__":
  args = sys.argv[1:]
  del sys.argv[1:]
  run(*tuple(args))
