import os, sys
import cProfile
from iotbx import pdb

def run(args):
  print "run the following arge: ",args
  print type(args)
  pdb_inp = pdb.input(args[0])
  pdb_obj = pdb.hierarchy.input(file_name=args[0])
  hierarchy = pdb_inp.construct_hierarchy()
  inpparameters = dir(pdb_inp)[22:]
  hparameters = dir(hierarchy)[23:]
  help(pdb_inp)
  print '|||'*50
  help(hierarchy)
  print '|||'*50
  help(pdb_obj)
  BIOMT_info = pdb_inp.process_BIOMT_records()
  #for i in range(0,len(inpparameters),4):
    #print inpparameters[i:i+4]
  #print '='*60
  #for i in range(0,len(hparameters),4):
    #print hparameters[i:i+4]  
  
  
  #print inpparameters[i:]
  #print inpparameters[:10]
  #print len(inpparameters)
  #print dir(hierarchy)
  print 'End run'

if __name__=="__main__":
  #args = sys.argv[1:]	# Use this line when working from the command line
  #del sys.argv[1:]
  args = ['1JS9.pdb'] 
  #cProfile.run('run(*tuple(args))')
  run(tuple(args))

  

