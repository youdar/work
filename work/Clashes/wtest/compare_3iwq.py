from __future__ import division
import os,sys


def print_rec_probe(lst):
    lst = lst.split(',')
    key1 = lst[0][2:-1]
    key2 = lst[1][2:-1]
    seq_i = lst[4]
    seq_j = lst[5]
    model = float(lst[6])
    vdw = float(lst[7])
    delta = model -vdw
    outstr = '{0:16} {1:16} {2:7}  {3:7} {4:6.2f} {5:6.2f} {6:6.2f}'.format(key1,key2, seq_i,seq_j,model,vdw,delta)
    print outstr





if __name__=='__main__':
    
    directory_path = 'c:\Phenix\Dev\Work\work\junk'
    directory_path = os.path.realpath(directory_path)
    internal_file_name = '3iwq-clashlist.txt'
    probe_file_name = '3iwq-prob.txt'
    probe = open(os.path.join(directory_path, probe_file_name), "r").readlines()
    internal = open(os.path.join(directory_path, internal_file_name), "r").readlines()
    
    for rec in probe:
        print_rec_probe(rec)
        
    #for rec in internal:
        #print_rec_probe(rec)
    
    


