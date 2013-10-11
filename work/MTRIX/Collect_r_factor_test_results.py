from __future__ import division
import os

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
        if (len(d) > 1):
            msg = '  '.join(d[-2:])
            msg = msg.replace('\n','  ',10)
            probelm_files.append(file + ':100:' + msg)                   
        elif (d == []):
            #print 'problem with the file {}'.format(file)
            probelm_files.append(file + ':100:')
        else:
            # remove '\n' from the end of the file name
            d = d[0].strip()
            if d.startswith('Sorry'):
                d = d.replace('\n','  ',10)
                probelm_files.append(file + ":100:" + d)
            else:
                t = d.split(':')
                pdb_file = t[0]
                r = float(t[1])                
                msg = ''.join(t[2:])
                #[pdb_file,r,msg] = d.split('::')
                if (msg == '') and (r < 10):
                    data.append([pdb_file,r,msg])
                    data_files.append('{0}::{1}::{2}'.format(file,t[1],'OK'))
                else:
                    probelm_files.append('{0}::{1}::{2}'.format(file,t[1],msg))
    print 'number of files with one data line: {}'.format(len(data))
    print 'number of files with problems: {}'.format(len(probelm_files))
    return data_files, probelm_files

def add_to_data_files(data_files, probelm_files,write_files=False):
    if write_files:
        f = open('Collect_tested_files',"a")
        g = open('files_with_problems',"a") 
        # write the results of good files
        for d in data_files:
            f.write(d + '\n')
        # write the results of files with issues
        for d in probelm_files:
            g.write(d + '\n')
        f.close()
        g.close()      
        
        
        
if __name__=='__main__':
    # locate the directory containing the log files
    directory_path = 'c:\Phenix\Dev\Work\work\queue_job'
    # convert the path to python format
    directory_path = os.path.realpath(directory_path)
    data_files, probelm_files = run(directory_path)
    print os.getcwd()
    #os.chdir('c:\\Phenix\\Dev\\Work\\work')
    #add_to_data_files(data_files, probelm_files,write_files=True)
    