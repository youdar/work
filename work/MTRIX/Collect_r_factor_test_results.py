from __future__ import division
#import numpy as np
#import matplotlib.mlab as mlab
#import matplotlib.pyplot as plt
import pylab as plb
import cPickle as pickle
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
        d = open(os.path.join(directory_path, file), "r").read().splitlines()
        file = file[4:]
        # check every line in d if it can be splited
        msg = []
        l_parts = d[0].split('::')
        file_ok = False
        if len(d) == 1:
            # Check if the log file starts with the file name
            if l_parts[0] == file:
                # check if this is the line with the r values
                # check if score is OK
                if l_parts[2]<'100':
                    # R value < 100 means that the file was processed correctly
                    temp_data_files = [file,float(l_parts[1]),float(l_parts[2]),float(l_parts[3])]
                    file_ok = True
                else:
                    # file with issues
                    file_ok = False

            # add information on rotation matix problems
            elif l_parts[0].startswith('Rotation matrices are not proper!'):
                msg.append('::Rotation matrices are not proper!')
                file_ok = False
            else:
                msg.append('::File was not processed properly!')
                file_ok = False
        # add errors from first line in file
        msg.append('::'.join(l_parts[4:]))
        if file_ok:
            temp_data_files.extend(msg)
            data_files.append(temp_data_files)
        else:
            temp_data_files = [file,msg]
            probelm_files.append(temp_data_files)

    print 'total number of files to process: {}'.format(len(files))
    print 'number of files with good data line: {}'.format(len(data_files))
    print 'number of files with problems: {}'.format(len(probelm_files))
    return data_files, probelm_files

def add_to_data_files(data_files, probelm_files,write_files=False):
    osType = sys.platform
    if osType.startswith('win'):
        directory_path = 'c:\Phenix\Dev\Work\work\MTRIX\Data'
    else:
        directory_path = '/net/cci-filer2/raid1/home/youval/Work/work/MTRIX/Data'
    os.chdir(directory_path)
    print 'add data in {}'.format(os.getcwd())
    if write_files:
        if os.path.isfile('Collect_tested_files'):
            Collect_tested_files = pickle.load(open('Collect_tested_files',"r"))
            Collect_tested_files.extend(data_files)
            pickle.dump(Collect_tested_files, open('Collect_tested_files','w'))
        else:
            pickle.dump(data_files, open('Collect_tested_files','w'))
        pickle.dump(probelm_files, open('files_with_problems','w'))



if __name__=='__main__':
    # locate the directory containing the log files
    osType = sys.platform
    if osType.startswith('win'):
        directory_path = r'c:\Phenix\Dev\Work\work\junk\queue_job_3'
    else:
        directory_path = '/net/cci-filer2/raid1/home/youval/Work/work/junk/queue_job_3'
    print os.getcwd()
    # convert the path to python format
    directory_path = os.path.realpath(directory_path)
    data_files, probelm_files = run(directory_path)
    #
    #add_to_data_files(data_files, probelm_files,write_files=True)
