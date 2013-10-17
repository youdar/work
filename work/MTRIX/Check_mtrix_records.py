from __future__ import division
import os,sys
from iotbx import pdb
import cPickle as pickle
from scitbx import matrix



def run():
    osType = sys.platform
    if osType.startswith('win'):
        directory_path = 'c:\Phenix\Dev\Work\work'
    else:
        directory_path = '/net/cci-filer2/raid1/home/youval/Work/work'
    print os.getcwd()
    os.chdir(directory_path)
    f = open('mtrix_not_included_in_pdb.txt','w')
    good_MTRIX_pdb_files = pickle.load(open('dict_good_MTRIX_pdb_files','r'))
    MTRIX_with_Straucture_Factor = pickle.load(open('MTRIX_with_Straucture_Factor_file_list','r'))
    results = []
    n = len(MTRIX_with_Straucture_Factor)
    for file_name in MTRIX_with_Straucture_Factor:
        pdb_inp = pdb.input(file_name=good_MTRIX_pdb_files[file_name])          # read the pdb file data
        TRASFORM_info = pdb_inp.process_mtrix_records(error_handle=False,eps=1e-2)
        TRASFORM_info = [matrix.rt(x.values) for x in TRASFORM_info if not x.coordinates_present]
        if TRASFORM_info != []:
            f.write(file_name + '\n')
            print '*'*60
            results.append(file_name)
    f.close()
    print 'files with MTRIX records that are not already include in pdb file: {}'.format(len(results))
    print 'Total files with mtrix records: {}'.format(n)


if __name__=='__main__':
    run()