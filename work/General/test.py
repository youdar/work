from __future__ import division
from libtbx.utils import Sorry
from iotbx import pdb
from scitbx import matrix
import os, sys
import string
import math
import cProfile
from  iotbx.pdb.multimer_reconstruction import *







#def run(args):
  #if (len(args) == 0):
  #raise RuntimeError("Please specify one or more pdb file names.")
  #for file_name in args:
  #pdb_obj = iotbx.pdb.hierarchy.input(file_name=file_name)
  #pdb_obj.hierarchy.overall_counts().show()
  #x = iotbx.pdb.input(file_name=file_name)
  #xh = x.construct_hierarchy()
  #biomt = x.process_BIOMT_records()


#def run(args):
  #if (len(args) == 0):
  #raise RuntimeError("Please specify one or more pdb file names.")
  #for file_name in args:
  #pdb_obj = iotbx.pdb.hierarchy.input(file_name=file_name)
  #pdb_obj.hierarchy.overall_counts().show()
  #for model in pdb_obj.hierarchy.models():
  #for chain in model.chains():
  #for rg in chain.residue_groups():
    #print 'resid: "%s"' % rg.resid()
    #for ag in rg.atom_groups():
    #print '  altloc: "%s", resname: "%s"' % (ag.altloc, ag.resname)
    #for atom in ag.atoms():
    #print '    ', atom.name



  #print '~'*50

def chains_names(BIOMT_transform_numbers,nChains, unique_chain_names):
  #def _chains_names(self):
  ''' (int, int, set) -> dictionary

  Create a dictionary
  keys: a string made of chain_name + str(BIOMT_transform_number)
  values: two letters and digits string combination

  The total number of new chains can be large and chain names might repeat themselves several times in a pdb file.
  We want to be able to asign the same name chains with the same name that go the same BIOMT transformation

  Arguments:
  BIOMT_transform_numbers -- an integer. the number of BIOMT transformation operations in the pdb object
  nChains -- an integer. the number of chains as interpreted by pdb.hierarchy.input
  unique_chain_names -- a set. a set of unique chain names

  Returns:
  new_names -- a dictionary. {'A1': 'aa', 'A2': 'gq',....} map a chain name and
  a BIOMT transform number to a new chain name
  '''
  # creat list of character from which to assemble the list of names
  total_chains_number = BIOMT_transform_numbers*len(unique_chain_names)
  chr_number = int(math.sqrt(total_chains_number)) + 1		# the number of charater needed to produce new names
  chr_list = list(string.ascii_letters) + list(string.digits)	# build character list
  chr_list = chr_list[:chr_number]				# take only as many characters as needed
  dictionary_values = set([ x+y for x in chr_list for y in chr_list])
  dictinary_key = set([x+str(y) for x in unique_chain_names for y in range(1,BIOMT_transform_numbers+1)])

  new_names_dictionary  = {x:y for (x,y) in zip(dictinary_key,dictionary_values)}
  return new_names_dictionary

def write_to_file(f,x):
  f.write(x)

def run_test():
  os.chdir('/net/cci/youval/Work/')
  m = multimer('1JS9.pdb','ba')
  #x = m.sites_cart()
  #y = m.get_xyz2()
  #print x == y

if (__name__ == "__main__"):
  #run(tuple(['1JS9.pdb']))
  #cProfile.run("chains_names(6,4,{'A','B'})")
  #print chains_names(3,4,{'A','B'})

  #pdb_input_file_name = '1S58.pdb'
  #pdb_inp = pdb.input(file_name=pdb_input_file_name)
  #pdb_obj = pdb.hierarchy.input(file_name=pdb_input_file_name)
  #pdb_obj_new = pdb_obj.hierarchy.deep_copy()
  #crystal_symmetry=pdb_inp.crystal_symmetry
  #os.chdir('/net/cci-filer2/raid1/home/youval/Work/work')
  #f = open('zz_out.txt','w')
  #write_to_file(f,'hi\n')
  #write_to_file(f,'end')
  #f.close()

  cProfile.run('run_test()')
  print 'ok'
