from __future__ import division
import pylab as plb
import cPickle as pickle
import sys,os


def run(directory_path):
    '''
    Data collected from run of submit_Clashes_to_queue.py

    This will read the files from c:\Phenix\Dev\Work\work\Clashes\queue_clash
    and orgenized it as follow:

    data = list([score_with_hydrogen,score_without_hydrogen]...)
    data_files = list([score_with_hydrogen,score_without_hydrogen,experment_type,file_name]...)
    data_dict[file_name] = [score_with_hydrogen,score_without_hydrogen,experment_type]

    '''
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
            # for l in d:
            # pdb_file_name::score_with_hydrogen::score_without_hydrogen::experment_type
            raw_data = d[0].split('::')
            file_name = raw_data[0]
            score_with_hydrogen = float(raw_data[1])
            score_without_hydrogen = float(raw_data[2])
            experment_type = {x.strip() for x in raw_data[3:]}
            experment_type = ','.join(experment_type)

            data.append([score_with_hydrogen,score_without_hydrogen])
            data_files.append([score_with_hydrogen,score_without_hydrogen,experment_type,file_name])
            data_dict[file_name] = [score_with_hydrogen,score_without_hydrogen,experment_type]
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