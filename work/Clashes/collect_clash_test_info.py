from __future__ import division
import pylab as plb
import cPickle as pickle
import sys,os


def run(directory_path):
    # Read files in directory_path
    files = os.listdir(directory_path)
    # collect only the files that starts with log_
    files = [x for x in files if x.startswith('log_')]
    data = []
    data_files = []
    data_dict = {}
    for file in files:
        d = open(os.path.join(directory_path, file), "r").readlines()
        file = file[4:]
        # check every line in d if it can be splited
        if len(d)==1:
            for l in d:
                [file_name,clash_score] = l.split('::')
                clash_score = float(clash_score.strip())
                data.append(clash_score)
                data_files.append([clash_score,file_name])
                data_dict[file_name] = clash_score
    return data,data_files,data_dict


def add_to_data_files(data,data_files,data_dict,work_path):
    pickle.dump(data, open('pdb_clash_scores','w'))
    pickle.dump(data_files, open('pdb_clash_score_and_name','w'))
    pickle.dump(data_dict, open('pdb_clash_score_dict','w'))

if __name__=='__main__':
    # locate the directory containing the log files
    osType = sys.platform
    if osType.startswith('win'):
        directory_path = 'c:\Phenix\Dev\Work\work\Clashes\queue_clash'
        work_path = 'c:\Phenix\Dev\Work\work\Clashes'
    else:
        directory_path = '/net/cci/youval/Work/work/Clashes/queue_clash'
        work_path = '/net/cci/youval/Work/work/Clashes'
    # convert the path to python format
    directory_path = os.path.realpath(directory_path)
    os.chdir(work_path)
    data,data_files,data_dict = run(directory_path)
    #add_to_data_files(data,data_files,data_dict,work_path)