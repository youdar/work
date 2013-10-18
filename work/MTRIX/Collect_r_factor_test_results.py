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
        file = file[4:]
        # files with errors have more than one line
        ld =  len(d)
        # check every line in d if it can be splited
        msg = []
        for l in d:
            l_parts = l.split('::')
            if l_parts[0] == file:
                # check if this is the line with the r values
                if len(l_parts)>3:
                    l_parts = [x.strip() for x in l_parts]
                    # check if score is OK
                    if l_parts[2]<'100':
                        temp_data_files = [file,float(l_parts[1]),float(l_parts[2]),float(l_parts[3])]
                        msg.append('::'.join(l_parts[4:]))
                        #temp_data.append(float(min(l_parts[1:3])))
                    else:
                        # collect files with issues
                        new_record = [file,'::'.join(l_parts[4:])]
                        probelm_files.append(new_record)
            # add information on rotation matix problems
            elif l_parts[0].startswith('Rotation matrices are not proper!'):
                msg.append('::Rotation matrices are not proper!')
        temp_data_files.extend(msg)
        data_files.append(temp_data_files)


    print 'number of files with good data line: {}'.format(len(data_files))
    print 'number of files with problems: {}'.format(len(probelm_files))
    # plot results
    #xend = len(data)+1
    #x = range(1, xend)
    #plb.plot(x,data,'o')
    #plb.show()

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
        pickle.dump(data_files, open('Collect_tested_files','w'))
        pickle.dump(probelm_files, open('files_with_problems','w'))



if __name__=='__main__':
    # locate the directory containing the log files
    osType = sys.platform
    if osType.startswith('win'):
        directory_path = 'c:\Phenix\Dev\Work\work\queue_job'
    else:
        directory_path = '/net/cci-filer2/raid1/home/youval/Work/work/queue_job'
    print os.getcwd()
    # convert the path to python format
    directory_path = os.path.realpath(directory_path)
    data_files, probelm_files = run(directory_path)
    print os.getcwd()
    #os.chdir('c:\\Phenix\\Dev\\Work\\work')
    #add_to_data_files(data_files, probelm_files,write_files=True)
