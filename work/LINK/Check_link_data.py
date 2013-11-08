from __future__ import division
from iotbx import pdb
from mmtbx.chemical_components import get_type
from iotbx.pdb import common_residue_names_get_class
import cPickle as pickle
import numpy as nm
import os,sys

'''
Exploring the contant of the LINK records in PDB files

LINK Record Format
http://www.wwpdb.org/documentation/format33/sect6.html#LINK

COLUMNS         DATA TYPE      FIELD           DEFINITION
------------------------------------------------------------------------------------
 1 -  6         Record name    "LINK  "
13 - 16         Atom           name1           Atom name.
17              Character      altLoc1         Alternate location indicator.
18 - 20         Residue name   resName1        Residue  name.
22              Character      chainID1        Chain identifier.
23 - 26         Integer        resSeq1         Residue sequence number.
27              AChar          iCode1          Insertion code.
43 - 46         Atom           name2           Atom name.
47              Character      altLoc2         Alternate location indicator.
48 - 50         Residue name   resName2        Residue name.
52              Character      chainID2        Chain identifier.
53 - 56         Integer        resSeq2         Residue sequence number.
57              AChar          iCode2          Insertion code.
60 - 65         SymOP          sym1            Symmetry operator atom 1.
67 - 72         SymOP          sym2            Symmetry operator atom 2.
74 � 78         Real(5.2)      Length          Link distance

@author: Youval Dar
'''

########################################################################
class Link_records(object):
    """
    Link records list from a PDB file
    """
    def __init__(self):
        self.data = []
        self.poslist = [[0,6],[12,16],[16,17],[17,20],[21,22],[22,26],[26,27],[42,46],
                        [46,47],[47,50],[51,52],[52,56],[56,57],[59,65],[66,72],[74,78]]
        self.posName = ['LINK','name1','altLoc1','resName1','chainID1','resSeq1','iCode1',
                        'name2','altLoc2','resName2','chainID2','resSeq2','iCode2','sym1','sym2','Length','fileNmae']
        self.name_location = {'resName1':3,'resName2':9,'Length':15}



def records_process(record_items):
    '''(list of list of strings) -> sets

    check the unique data options for each column
    in the LINK record
    '''
    record_length = len(record_items.data[0]) - 1
    print 'There are {0} Records with {1} components'.format(len(record_items.data),(record_length-1))
    print '*'*50
    # the features we are interested in
    features_list = ['name1','altLoc1','resName1','chainID1','iCode1',
        'name2','altLoc2','resName2','chainID2','iCode2']
    # Collect the set of unique features
    for i in range(1,record_length):
        r = set(x[i] for x in record_items.data if record_items.posName[i] in features_list)
        if record_items.posName[i] in features_list:
            print 'elements set for {} is: '.format(record_items.posName[i])
            print '-'*50
            print_set(r)
            print '-'*50

def process_links(record_items):
    # Collect link types dictinaries
    link_file_dict = {}		# dictionary to files containing the LINK type
    link_type_dict = {}		# dictionary to the number of LINKs of particular type
    link_length_dict = {}	# dictionary to the length of LINKs of particular type
    # Open dictionary with number of models for each file
    file_name_to_model_number_dict = pickle.load(open('file_name_to_model_number_dict','r'))
    # Name locations
    n1 = record_items.name_location['resName1']
    n2 = record_items.name_location['resName2']
    n3 = record_items.name_location['Length']
    for x in record_items.data:
        # Collect resName1 and resName2 from each LINK record
        names = [x[n1],x[n2]]
        # get the name of the file containing the current link
        file_name = x[-1]
        link_length = x[n3]
        # Sort to ignor order
        names.sort()
        Res_type = map(get_res_type,names)
        # Process cases where get_type returns None
        Res_type = [Res_type[0].strip(),Res_type[1].strip()]		# remove spaces
        Res_type = [Res_type[0].strip('"'),Res_type[1].strip('"')]	# remove "
        Res_type = [Res_type[0].lower(),Res_type[1].lower()]		# convert all to lower case
        link_type = '{0:<33} - {1:<33}'.format(*Res_type)			# compose link
        if link_type_dict.has_key(link_type):
            link_type_dict[link_type] += 1
            link_file_dict[link_type] += [file_name]
            if link_length != ' ':
                link_length_dict[link_type] += [float(link_length)]
        else:
            link_type_dict[link_type] = 1
            link_file_dict[link_type] = [file_name]
            link_length_dict[link_type] = [float(link_length)]
    # print LINK type results
    print ' '*35 + 'Link type list'
    print ' '*35 + 'Number of unique links {}'.format(len(link_type_dict))
    print '{0:<69}{1:<9}{2:>6}{3:>9}'.format('Link Type','LINKS Number','  Length','File')
    print '_'*99
    for link_type,count in sorted(link_type_dict.items()):
        # look for first file with 1 model
        link_file_dict[link_type][0]
        n_models = 1000
        for rec in link_file_dict[link_type]:
            key = rec[3:7]	# take only the 4 letter pdb file name
            if key in file_name_to_model_number_dict:
                n_models_tmp = file_name_to_model_number_dict[key]
            else:
                n_models_tmp = 1000
            if n_models_tmp == 1:
                first_file = rec
                break
            elif (n_models_tmp > 0) and (n_models > n_models_tmp):
                n_models = n_models_tmp
                first_file = rec

        mean = sum(link_length_dict[link_type])/len(link_length_dict[link_type])
        print '{0:<76}{1:<9}{2:<6.2f}{3:>15}'.format(link_type,count,mean,first_file)
    # Save results dictionaries
    #pickle.dump(link_file_dict,open('LINK_type_to_file_dict','w'))
    #pickle.dump(link_type_dict,open('LINK_type_to_num_dict','w'))
    #pickle.dump(link_length_dict,open('LINK_type_to_length','w'))
    #pickle.dump(sorted(link_type_dict.iterkeys()),open('LINK_type_keys','w'))
    # To retrive the pickled data use
    # data = pickle.load(open('file_name','r'))

def get_res_type(res_name):
    '''(string) -> string
    process residue type by the residue name
    '''
    res_type = common_residue_names_get_class(res_name)
    if res_type == 'other':
        res_type = get_type(res_name)
        if res_type == None:
            res_type = 'None'
    return res_type


def print_set(s):
    '''
    Print s in a table with 8 columns
    '''
    s = list(s)
    l = len(s)
    # pad s with empty string to make the length devisible by 8
    s.extend(['']*(8 - l%8))
    for i in range(l//8+1):
        print '{0:<12}{1:<12}{2:<12}{3:<12}{4:<12}{5:<12}{6:<12}{7:<12} '.format(*s[i*8:(i+1)*8])
        # the * in the format is for unpacking the list

def run(raw_link_data):
    '''(list of strings) -> print info to screen


    poslist is the location of the records in LINK accroding to
    http://www.bmsc.washington.edu/CrystaLinks/man/pdb/part_48.html  (see remarks at the beginning
    of this file
    '''
    data_split = Link_records()
    num_of_links_in_file = []
    data_split.data = []
    record_counter = 0
    for line in raw_link_data:
        if line.startswith('LINK'):
            # get the file name
            l = line.split()[-1]
            # Split the LINK record
            if len(line) >= 93:
                b = [line[x[0]:x[1]] for x in data_split.poslist]
            else:
                b = [line[x[0]:x[1]] for x in data_split.poslist[:11]]
                b.extend([' ']*5)
            b.append(l.strip())			# add the file name
            record_counter += 1
            data_split.data.append(b)			# keep the split data
        else:
            num_of_links_in_file.append(record_counter)
            record_counter = 0
    print 'number of files with LINK records: {}'.format(len(num_of_links_in_file))
    print 'average number of records: {:.1f}'.format(nm.mean(num_of_links_in_file))
    print 'standard deviation of number of records: {:.1f}'.format(nm.std(num_of_links_in_file))
    print 'max number of records: {}'.format(nm.max(num_of_links_in_file))
    print 'min number of records: {}'.format(nm.min(num_of_links_in_file))
    print '='*50
    process_links(data_split)
    print '='*50
    records_process(data_split)
    print '='*50
    print 'done'


if __name__=='__main__':
    # work in my working directory
    osType = sys.platform
    if osType.startswith('win'):
        os.chdir('c:\Phenix\Dev\Work\work\LINK')
    else:
        os.chdir('/net/cci-filer2/raid1/home/youval/Work/work/LINK')
    # open the raw data file that was produced by gather_link_info.py
    f = open('Link_pdb_file_raw_data.txt','r')
    raw_link_data = f.readlines()
    f.close()
    run(raw_link_data)