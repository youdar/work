from __future__ import division
import pylab as plb
import cPickle as pickle
import sys,os


def run(directory_path):
    '''
    Data collected from run of submit_clashscore_Clashes_to_queue.py

    This will read the files from c:\Phenix\Dev\Work\work\Clashes\queue_clash_compare_...
    and orgenized it as follow:

    data = list([total_nb_clashscore,without_sym_nb_clashscore,clashscore_probe]...)
    data_files = list([total_nb_clashscore,without_sym_nb_clashscore,clashscore_probee,file_name]...)
    data_dict[file_name] = [total_nb_clashscore,without_sym_nb_clashscore,clashscore_probe]

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
        d = open(os.path.join(directory_path, file), "r").read().splitlines()
        file = file[4:]
        # check every line in d if it can be splited
        if len(d)==1 and d[0][4:6]=='::':
            # for l in d:
            # pdb_file_name::total_nb_clashscore::without_sym_nb_clashscore::clashscore_probe::run_time::run_time_probe
            raw_data = d[0].split('::')
            file_name = raw_data[0]
            total_nb_clashscore = float(raw_data[1])
            without_sym_nb_clashscore = float(raw_data[2])
            clashscore_probe = float(raw_data[3])
            total_nb_clashscore_time = float(raw_data[4])
            clashscore_probe_time = float(raw_data[5])
            #
            data.append([total_nb_clashscore,without_sym_nb_clashscore,clashscore_probe])
            data_files.append([total_nb_clashscore,without_sym_nb_clashscore,clashscore_probe,total_nb_clashscore_time,clashscore_probe_time,file_name])
            data_dict[file_name] = [total_nb_clashscore,without_sym_nb_clashscore,clashscore_probe]
        elif d != []:
            if len(d) > 1:
                #print file
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


def add_to_data_files(data,data_files,data_dict,out_directory,outFileName):
    filename = os.path.join(out_directory,outFileName)
    pickle.dump(data, open(filename,'w'))
    pickle.dump(data_files, open(filename+'_and_name','w'))
    pickle.dump(data_dict, open(filename+'_dict','w'))

def write_to_txt_file(data_files,out_directory,outFileName):
    '''
    lines in text file
    file_name,total_nb_clashscore,without_sym_nb_clashscore,clashscore_probe,total_nb_clashscore_time,clashscore_probe_time
    '''
    filename = os.path.join(out_directory,outFileName)
    f = open(filename+'.txt','w')
    for x in data_files:
        outstr = '{5},{0:.2f},{1:.2f},{2:.2f},{3:.1f},{4:.1f}\n'.format(*x)
        f.write(outstr)
    f.close()

if __name__=='__main__':
    # locate the directory containing the log files
    queue_folder = 'queue_clash_compare_12_11_2013'
    outFileName = 'clashscore_compare_reduce_12_11_2013'
    outflder = 'Data'
    osType = sys.platform
    if osType.startswith('win'):
        directory_path = 'c:\Phenix\Dev\Work\work\Clashes'
        work_path = 'c:\Phenix\Dev\Work\work\Clashes'
    else:
        directory_path = '/net/cci/youval/Work/work/Clashes'
        work_path = '/net/cci/youval/Work/work/Clashes'
    directory_path = os.path.join(directory_path,queue_folder)
    out_directory = os.path.join(work_path,outflder)
    # convert the path to python format
    directory_path = os.path.realpath(directory_path)
    #os.chdir(work_path)
    data,data_files,data_dict,files_with_problems = run(directory_path)
    add_to_data_files(data,data_files,data_dict,out_directory,outFileName)
    write_to_txt_file(data_files,out_directory,outFileName)
    print len(data_files)
    print 'Done...'