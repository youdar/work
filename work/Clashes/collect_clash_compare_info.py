from __future__ import division
import pylab as plb
import cPickle as pickle
import sys,os


def run(directory_path):
    '''
    Data collected from run of submit_clashscore_Clashes_to_queue.py

    This will read the files from c:\Phenix\Dev\Work\work\Clashes\queue_clash_compare
    and orgenized it as follow:

    data = list([clashscore_internal,lashscore_probe]...)
    data_files = list([clashscore_internal,lashscore_probe,file_name]...)
    data_dict[file_name] = [sclashscore_internal,lashscore_probe]

    '''
    # Read files in directory_path
    files = os.listdir(directory_path)
    # collect only the files that starts with log_
    files = [x for x in files if x.startswith('log_')]
    data = []
    data_files = []
    data_dict = {}
    files_with_problems = []
    test = set()
    for file in files:
        d = open(os.path.join(directory_path, file), "r").readlines()
        file = file[4:]
        # check every line in d if it can be splited
        if len(d)==1 and d[0][4:6]=='::':
            # for l in d:
            # pdb_file_name::clashscore_internal::lashscore_probe,file_name
            raw_data = d[0].split('::')
            file_name = raw_data[0]
            clashscore_internal = float(raw_data[1])
            clashscore_probe = float(raw_data[2])
            clashscore_internal_time = float(raw_data[3])
            clashscore_probe_time = float(raw_data[4])
            #
            data.append([clashscore_internal,clashscore_probe])
            data_files.append([clashscore_internal,clashscore_probe,clashscore_internal_time,clashscore_probe_time,file_name])
            data_dict[file_name] = [clashscore_internal,clashscore_probe]
        elif d != []:
            if len(d) > 1:
                print file
                for x in d:
                    if 'KeyError' in x:
                        new_key = x.split(':')[-1].strip()
                        test.add(new_key)
                        files_with_problems.append([file,new_key])
            else:
                files_with_problems.append([file,d[0]])

        else:
            files_with_problems.append([file_name,[]])
    return data,data_files,data_dict,files_with_problems


def add_to_data_files(data,data_files,data_dict,work_path):
    pickle.dump(data, open('clashscore_compare_ref_1','w'))
    pickle.dump(data_files, open('clashscore_compare_and_name_ref_1','w'))
    pickle.dump(data_dict, open('clashscore_compare_dict_ref_1','w'))

def write_to_txt_file(data_files):
    f = open('clashscore_compare_and_name_ref_1.txt','w')
    for x in data_files:
        outstr = '{4},{0:.1f},{1:.1f},{2:.1f},{3:.1f}\n'.format(*x)
        f.write(outstr)
    f.close()

if __name__=='__main__':
    # locate the directory containing the log files
    osType = sys.platform
    if osType.startswith('win'):
        directory_path = 'c:\Phenix\Dev\Work\work\Clashes\queue_clash_compare'
        work_path = 'c:\Phenix\Dev\Work\work\Clashes'
    else:
        directory_path = '/net/cci/youval/Work/work/Clashes/queue_clash_compare'
        work_path = '/net/cci/youval/Work/work/Clashes'
    # convert the path to python format
    directory_path = os.path.realpath(directory_path)
    os.chdir(work_path)
    data,data_files,data_dict,files_with_problems = run(directory_path)
    #add_to_data_files(data,data_files,data_dict,work_path)
    #write_to_txt_file(data_files)
    print len(data_files)
    print 'Done...'