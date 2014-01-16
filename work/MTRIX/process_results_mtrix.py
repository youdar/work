from __future__ import division
#import numpy as np
#import matplotlib.mlab as mlab
#import matplotlib.pyplot as plt
import cPickle as pickle
import pylab as plb
import os, sys


def run(directory_path):
    # Read files in directory_path
    files = os.listdir(directory_path)
    # collect only the files that starts with log_
    files = [x for x in files if x.startswith('log_')]
    data = []
    data_files = []
    probelm_files = []
    for file in files:
        d = open(os.path.join(directory_path, file), "r").readlines()
        file = file[4:]
        if d == []:
            data_files.append(file)
        else:
            probelm_files.append(file)
    print len(probelm_files)
    print len(data_files)

    return data_files, probelm_files

def add_to_data_files(data_files, probelm_files,write_files=False):
    osType = sys.platform
    if osType.startswith('win'):
        directory_path = 'c:\Phenix\Dev\Work\work'
    else:
        directory_path = '/net/cci-filer2/raid1/home/youval/Work/work'
    os.chdir(directory_path)
    print os.getcwd()
    if write_files:
        pickle.dump(data_files, open('files_with_good_MTRIX','w'))
        pickle.dump(probelm_files, open('files_with_bad_MTRIX','w'))

if __name__=='__main__':
    # locate the directory containing the log files
    osType = sys.platform
    if osType.startswith('win'):
        directory_path = 'c:\Phenix\Dev\Work\work\queue_job_2'
    else:
        directory_path = '/net/cci-filer2/raid1/home/youval/Work/work/queue_job_2'
    print os.getcwd()
    # convert the path to python format
    directory_path = os.path.realpath(directory_path)
    data_files, probelm_files = run(directory_path)
    print os.getcwd()
    #os.chdir('c:\\Phenix\\Dev\\Work\\work')
    #add_to_data_files(data_files, probelm_files,write_files=True)