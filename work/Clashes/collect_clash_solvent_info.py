from __future__ import division
import pylab as plb
import cPickle as pickle
import sys,os


def run(directory_path):
    '''
    Data collected from run of submit_clashscore_Clashes_to_queue.py

    This will read the files from c:\Phenix\Dev\Work\work\Clashes\queue_clash_solvent_...
    and orgenized it as follow:

    data = list([clashscore_all_clashes,clashscore_simple,clashscore_only_sym_op,clashscore_solvent_solvent]...)
    data_files = list([clashscore_all_clashes,clashscore_simple,clashscore_only_sym_op,clashscore_solvent_solvent,file_name]...)
    data_dict[file_name] = [clashscore_all_clashes,clashscore_simple,clashscore_only_sym_op,clashscore_solvent_solvent]

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
            # pdb_file_name::clashscore_all_clashes::clashscore_simple::clashscore_only_sym_op::clashscore_solvent_solvent
            raw_data = d[0].split('::')
            file_name = raw_data[0]
            clashscore_all_clashes = float(raw_data[1])
            clashscore_simple = float(raw_data[2])
            clashscore_only_sym_op = float(raw_data[3])
            clashscore_solvent_solvent = float(raw_data[4])
            #
            data.append([clashscore_all_clashes,clashscore_simple,clashscore_only_sym_op,clashscore_solvent_solvent])
            data_files.append([clashscore_all_clashes,clashscore_simple,clashscore_only_sym_op,clashscore_solvent_solvent,file_name])
            data_dict[file_name] = [clashscore_all_clashes,clashscore_simple,clashscore_only_sym_op,clashscore_solvent_solvent]
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
    file_name,clashscore_all_clashes,clashscore_simple,clashscore_only_sym_op,clashscore_solvent_solvent
    '''
    filename = os.path.join(out_directory,outFileName)
    f = open(filename+'.txt','w')
    print filename+'.txt'
    for x in data_files:
        outstr = '{4},{0:.2f},{1:.2f},{2:.2f},{3:.1f}\n'.format(*x)
        f.write(outstr)
    f.close()

if __name__=='__main__':
    # locate the directory containing the log files
    queue_folder = 'queue_clash_solvent_12_13_2013'
    outFileName = 'clashscore_solvent_reduce_12_13_2013'
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
    #add_to_data_files(data,data_files,data_dict,out_directory,outFileName)
    #write_to_txt_file(data_files,out_directory,outFileName)
    print len(data_files)
    print 'Done...'