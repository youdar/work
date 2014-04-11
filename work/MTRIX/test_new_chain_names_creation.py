from __future__ import division
from pprint import pprint
import string
import math

__author__ = 'Youval'



def test_new_chain_names(i_transforms, unique_chain_names):
  ''' (int, int, set) -> dictionary

  Create a dictionary
  keys: a string made of chain_name + str(i_transform number)
  values: string (one or two chr long). letters and digits string combination

  The total number of new chains can be large (order of hundreds)
  Chain names might repeat themselves several times in a pdb file
  We want copies of chains with the same name to still have the same name after
  similar BIOMT/MTRIX transformation

  Arguments:
  i_transform : (integer) the number of BIOMT/MTRIX transformation operations
                in the pdb object
  unique_chain_names : (set) a set of unique chain names

  Returns:
  new_names : a dictionary. {'A1': 'aa', 'A2': 'gq',....} map a chain name and
  a transform number to a new chain name

  >>> self._chains_names(3,4,{'A','B'})
  {'A1': 'aa', 'A3': 'ac', 'A2': 'ab', 'B1': 'ba', 'B2': 'bb', 'B3': 'bc'}
  '''
  # create list of character from which to assemble the list of names
  total_chains_number = i_transforms*len(unique_chain_names)
  n_unique_chains = len(unique_chain_names)
  # start naming chains with a single letter
  chr_list1 = list(set(string.ascii_uppercase) - unique_chain_names)
  chr_list2 = list(set(string.ascii_lowercase) - unique_chain_names)
  chr_list1.sort()
  chr_list2.sort()
  new_names_list = chr_list1 + chr_list2
  # check if we need more chain names
  if len(new_names_list) < total_chains_number:
    n_names =  total_chains_number - len(new_names_list)
    # the number of character needed to produce new names
    chr_number = int(math.sqrt(n_names)) + 1
    # build character list
    chr_list = list(string.ascii_uppercase) + \
               list(string.ascii_lowercase) + \
               list(string.digits)
    # take only as many characters as needed
    chr_list = chr_list[:chr_number]
    extra_names = set([ x+y for x in chr_list for y in chr_list])
    # make sure not using existing names
    extra_names = list(extra_names - unique_chain_names)
    extra_names.sort()
    new_names_list.extend(extra_names)
  assert len(new_names_list) >= total_chains_number
  dictionary_values = new_names_list[:total_chains_number]
  dictinary_key = list(set([x+str(y+1) for x in unique_chain_names for y in
                       range(i_transforms)]))
  dictinary_key.sort()
  # create the dictionary
  zippedlists = zip(dictinary_key,dictionary_values)
  new_names_dictionary ={x:y for (x,y) in zippedlists}
  return new_names_dictionary

if __name__=='__main__':
  for (i_transforms, unique_chain_names) in [
    (0,{'A'}),(0,set()),(1,set()),(10,{'A','B'}),(10,{'A','B','D'}),
    (30,{'A','B','D'})]:
    outstr = 'i_transforms: {0} , unique_chain_names: {1}'
    print outstr.format(i_transforms, unique_chain_names)
    pprint(test_new_chain_names(i_transforms, unique_chain_names))
    print '-'*80