from __future__ import division
from iotbx import pdb
from mmtbx.chemical_components import get_type
from iotbx.pdb import common_residue_names_get_class
import cPickle as pickle
import numpy as nm
import os,sys


'''
go over the list of all files with LINK records and check how many models they have

record the result in a dictionary file_name_to_model_number_dict

If file does not exist, set the number of models to -1.
So when we look at the dictionary there will be some entry for it
'''

def run (files_list):
    print 'Start'
    print '*'*50
    file_name_to_model_number_dict = {}
    for file_name in files_list:
        file_name = file_name.strip()			# Remove '\n' from the end of the file name
        if os.path.isfile(file_name):
            pdb_inp = pdb.input(file_name=file_name)
            n_models = len(pdb_inp.model_ids())		# Get the number of models
            file_name = get_file_name(file_name)
            file_name_to_model_number_dict[file_name] = n_models
        else:
            file_name_to_model_number_dict[file_name] = -1
        print file_name
    pickle.dump(file_name_to_model_number_dict,open('file_name_to_model_number_dict','w'))
    print '*'*50
    print 'Done'

def get_file_name(file_name):
    '''(str)  -> str
    clean a pdb file name, remove path and file extensions

    Return the 4 letter pdb id
    '''
    osType = sys.platform
    if osType.startswith('win'):
        file_name = file_name.split('\\')[-1]
    else:
        file_name = file_name.split('/')[-1]
    file_name = file_name.split('.')[0]
    if len(file_name)>4:
        if 'pdb' in file_name:
            i = file_name.find('pdb')
            file_name = file_name[i+3:i+7]
    return file_name


if __name__=='__main__':
    osType = sys.platform
    if osType.startswith('win'):
        print 'Please run at LBL Chevy so that files can be handle localy'
    else:
        os.chdir('/net/cci-filer2/raid1/home/youval/Work/work/LINK')
    # open the raw data file that was produced by gather_link_info.py
    files_list = open('LINK_files_with_link_records.txt','r').readlines()
    print 'nuber of files to process: {}'.format(len(files_list))
    run(files_list)