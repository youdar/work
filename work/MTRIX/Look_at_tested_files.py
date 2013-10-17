from __future__ import division
import numpy as np
from numpy.random import rand
#import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import cPickle as pickle
import pylab as plb
import os, sys


def get_data():
    tested_files = open('Collect_tested_files',"r").read().splitlines()
    tested_files = [x.split('::') for x in tested_files]
    tested_files = [[x[0],float(x[1]),x[2]] for x in tested_files]
    r_factors_only = [x[1] for x in tested_files]
    return tested_files,r_factors_only


    
def plot_data(data,data_labels):
    
    def onpick3(event):
        ind = event.ind
        i = ind[0]
        print '{0}  R-Work: {1:.3f} (Model-Work)/Work: {2:.3f} '.format(data_labels[i],data[i],data[i])
        #print 'onpick3 scatter:', ind, npy.take(x, ind), npy.take(y, ind)  
    
    gr = 1.61803398875
    h = 10
    w = gr*h    
    d = 0.03
    fig = plt.figure(figsize=(w,h)) 
    plt.subplots_adjust(left=d, right=1-d, top=1-d, bottom=d)
    ax1 = fig.add_subplot(111)
    #col = ax1.scatter(z, y, 100*s, c, picker=True)
    l = len(data)
    x = xrange(1,l+1)
    s = 0.5 + abs(rand(l) - 1)
    c = rand(l)
    col = ax1.scatter(x, data, 100*s, c, picker=True)
    #fig.savefig('pscoll.eps')
    fig.canvas.mpl_connect('pick_event',onpick3)  
    fig.set_size_inches(w,h)
    ax1.set_ylim([0,max(data)+0.2])
    ax1.set_xlim([0,l+1])
    
    plt.show()    

def hist_data(data):
    plt.figure(2)
    plt.hist(data,bins=20)
    plt.show()    

def get_sets(data_files):
    all_files = {x[0] for x in data_files}
    equal_r = {x[0] for x in data_files if x[1]==x[2]} 
    return all_files,equal_r

def show_dara_summeries(data_files,all_files,equal_r,
                        not_included_mtrix,files_with_bad_MTRIX,
                        files_with_good_MTRIX):
    first_or_second = [1*(x[1]<x[2]) for x in data_files if x[1]!=x[2]]
    not_equal_r = all_files - equal_r
    r01 = {x[0] for x in data_files if min(x[1],[2]) <= 0.1}
    r02 = {x[0] for x in data_files if min(x[1],[2]) > 0.1 and min(x[1],[2]) <= 0.2 }
    r03 = {x[0] for x in data_files if min(x[1],[2]) > 0.2 and min(x[1],[2]) <= 0.3 }
    rrest = {x[0] for x in data_files if min(x[1],[2]) > 0.3}
    print 'Number of files that the recondtructes R-value is different than the pdb one'
    print len(not_equal_r)
    print 'Number of files that the reconstructed is smaller'
    print sum(first_or_second)
    print 'MTRIX transform are not included in pdb file data  '
    print '*'*60
    print '      value <= 0.1 {}'.format(len(not_included_mtrix & r01))
    print '0.1 < value <= 0.2 {}'.format(len(not_included_mtrix & r02))
    print '0.2 < value <= 0.3 {}'.format(len(not_included_mtrix & r03))
    print '      value >  0.3 {}'.format(len(not_included_mtrix & rrest))
    print 'Files with bad rotation MTRIX  '
    print '*'*60
    print '      value <= 0.1 {}'.format(len(files_with_bad_MTRIX & r01))
    print '0.1 < value <= 0.2 {}'.format(len(files_with_bad_MTRIX & r02))
    print '0.2 < value <= 0.3 {}'.format(len(files_with_bad_MTRIX & r03))
    print '      value >  0.3 {}'.format(len(files_with_bad_MTRIX & rrest))
    print 'Files with good rotation MTRIX  '
    print '*'*60
    print '      value <= 0.1 {}'.format(len(files_with_good_MTRIX & r01))
    print '0.1 < value <= 0.2 {}'.format(len(files_with_good_MTRIX & r02))
    print '0.2 < value <= 0.3 {}'.format(len(files_with_good_MTRIX & r03))
    print '      value >  0.3 {}'.format(len(files_with_good_MTRIX & rrest))
    print '*'*60
    a = not_equal_r & not_included_mtrix
    b = not_equal_r - a
    c = not_included_mtrix - a
    print len(a)
    print len(b)
    print len(c)
    print '*'*60    

if __name__=='__main__':
    # locate the directory containing the log files
    osType = sys.platform
    if osType.startswith('win'):
        directory_path = 'c:\Phenix\Dev\Work\work'
    else:
        directory_path = '/net/cci-filer2/raid1/home/youval/Work/work'
    # convert the path to python format
    directory_path = os.path.realpath(directory_path)
    os.chdir(directory_path)
    #print os.getcwd()
    data_files = pickle.load(open('Collect_tested_files','r'))
    probelm_files = pickle.load(open('files_with_problems','r'))    
    files_with_good_MTRIX = set(pickle.load(open('files_with_good_MTRIX','r')))
    files_with_bad_MTRIX = set(pickle.load(open('files_with_bad_MTRIX','r')))    
    not_included_mtrix = set(open('mtrix_not_included_in_pdb.txt', 'r').read().splitlines())
    # plot best data
    best_data = [min(x[1:]) for x in data_files]
    data_labels = [x[0] for x in data_files]
    # create sets
    all_files,equal_r = get_sets(data_files)
    show_dara_summeries(data_files,all_files,equal_r,
                        not_included_mtrix,files_with_bad_MTRIX,
                        files_with_good_MTRIX)
    plot_data(best_data, data_labels)
    #hist_data(best_data)    
    print 'done'
