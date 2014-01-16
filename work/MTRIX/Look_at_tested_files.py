from __future__ import division
import numpy as np
from numpy.random import rand
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
        print '{0}  R-Work: {1:.3f}'.format(data_labels[i],data[i])
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

def plot_data2(data):

    gr = 1.61803398875
    h = 10
    w = gr*h
    d = 0.03
    fig = plt.figure(figsize=(w,h))
    plt.subplots_adjust(left=d, right=1-d, top=1-d, bottom=d)
    ax1 = fig.add_subplot(111)
    #col = ax1.scatter(z, y, 100*s, c, picker=True)
    data.sort()
    l = len(data)
    x = xrange(1,l+1)
    s = 0.5
    c = rand(l)
    #col = ax1.scatter(x,data)
    col = ax1.scatter(x,data, 100*s, c)
    #fig.savefig('pscoll.eps')
    #fig.set_size_inches(w,h)
    #ax1.set_ylim([0,max(data)+0.2])
    #ax1.set_xlim([0,l+1])
    #plt.show()

    plt.figure(2)
    plt.hist(data,bins=100)
    plt.title('Compare R-work, (reported - calculated)/calculated')
    #plt.title('Compare R-work, (reconstructed - expected)/expected')
    plt.xlim([-1,1])
    plt.show()

def plot_data3(data):

    gr = 1.61803398875
    h = 10
    data.sort()
    x = [r[0] for r in data]
    y = [r[1] for r in data]
    l = len(data)
    p = plt.plot(x,y,'.')
    plt.xlabel('R-work ratio, (reconstructed - expected)/expected')
    plt.ylabel('Expected R-work')
    plt.show()

def plot_data4(data_NCS,data_ASU,data_issues,plt_lim=0.7):

    #gr = 1.61803398875
    gr = 1
    h = 10
    w = gr*h
    d = 0.08
    fig = plt.figure(figsize=(w,h))
    plt.subplots_adjust(left=d, right=1-d, top=1-d, bottom=d)
    ax1 = fig.add_subplot(111)
    #col = ax1.scatter(z, y, 100*s, c, picker=True)
    data_NCS.sort()
    data_ASU.sort()
    #l = len(data)
    #x = xrange(1,l+1)
    #s = [1000*rec[1] for rec in data]
    s = 40
    #c = rand(l)
    c_ncs = 'b'
    c_asu = 'y'
    c_issues = 'black'
    x_ncs = [rec[0] for rec in data_NCS]
    y_ncs = [rec[1] for rec in data_NCS]
    x_asu = [rec[0] for rec in data_ASU]
    y_asu = [rec[1] for rec in data_ASU]
    x_issues = [rec[0] for rec in data_issues]
    y_issues = [rec[1] for rec in data_issues]
    #col = ax1.scatter(x,data)
    p1 = ax1.scatter(x_asu,y_asu, s, c_asu, linewidths=0)
    p2 = ax1.scatter(x_ncs,y_ncs, s, c_ncs, linewidths=0)
    p3 = ax1.scatter(x_issues,y_issues, s, c_issues, linewidths=0)
    #fig.savefig('pscoll.eps')
    #fig.set_size_inches(w,h)
    #ax1.set_ylim([0,max(data)+0.2])
    #ax1.set_xlim([0,l+1])
    plt.xlim([0,plt_lim])
    plt.ylim([0,plt_lim])
    #
    #plt.xlabel('R-work calculated from PDB file',fontsize=18)
    #plt.ylabel('R-work calculated from reconstructed PDB file',fontsize=18)
    #plt.legend([p1,p2,p3],['From PDB','Reconstructed','R-work discrepancies'])
    #
    plt.xlabel('R-work calculated from PDB file - R-work reported',fontsize=18)
    plt.ylabel('R-work calculated from reconstructed PDB file - R-work reported',fontsize=18)
    plt.legend([p1,p2],['From PDB','Reconstructed'])
    #
    #plt.xlabel('R-work from PDB file REMARKS',fontsize=18)
    #plt.ylabel('Calculated R-work ',fontsize=18)
    #plt.legend([p1,p2],['Reported in PDB','Reconstructed'])
    #
    plt.show()

def hist_data(data):
    plt.figure(2)
    plt.hist(data,bins=20)
    plt.show()

def get_sets(data_files):
    all_files = {x[0] for x in data_files}
    equal_r = {x[0] for x in data_files if x[1]==x[2]}
    return all_files,equal_r

def show_data_summeries(data_files,all_files,equal_r,
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
    print '# files_with_good_MTRIX: {}'.format(len(files_with_good_MTRIX))
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
        directory_path = r'c:\Phenix\Dev\Work\work\MTRIX\Data'
    else:
        directory_path = '/net/cci-filer2/raid1/home/youval/Work/work/MTRIX/Data'
    # convert the path to python format
    directory_path = os.path.realpath(directory_path)
    os.chdir(directory_path)
    print os.getcwd()
    data_files = pickle.load(open('Collect_tested_files','r'))
    # data_files is a list of records [['3kk6', 0.207, 0.312, 0.312, 'OK'],...
    # ['3kk6', 0.207, 0.312, 0.312, 'OK'] : [pdb_file_name,r_work_expected,r_work_model_pdb,r_work_model_reconstructed,processing_msg]
    probelm_files = pickle.load(open('files_with_problems','r'))
    files_with_good_MTRIX = set(pickle.load(open('files_with_good_MTRIX','r')))
    files_with_bad_MTRIX = set(pickle.load(open('files_with_bad_MTRIX','r')))
    not_included_mtrix = set(open('mtrix_not_included_in_pdb.txt', 'r').read().splitlines())
    # plot best data
    #best_data = [min(x[1:]) for x in data_files]
    #data_labels = [x[0] for x in data_files]
    # create sets
    #all_files,equal_r = get_sets(data_files)
    #show_data_summeries(data_files,all_files,equal_r,
                        #not_included_mtrix,files_with_bad_MTRIX,
                        #files_with_good_MTRIX)
    # records in data_files
    # file_name::r_work_expected::r_work_model_pdb::r_work_model::msg
    #x[0]: file_name: 4 charaters PDB file name
    #x[1]: r_work_expected: r_work from pdb file
    #x[2]: r_work_model_pdb: r_work calulated from pdb
    #x[3]: r_work_model:r_work calculated for complete ASU
    #plot_data(best_data, data_labels)
    #hist_data(best_data)
    # plot (R_reconstructed-R_expected)/R_expected
    #data = [(x[1]-x[3])/x[3] for x in data_files]
    #data = [(x[3]-x[2]) for x in data_files ]
    #plot_data2(data)
    #
    # seperate between the files with complete ASU and a single NCS
    NCS = []
    ASU = []
    for x in data_files:
        if x[0] in not_included_mtrix:
            NCS.append(x)
        else:
            ASU.append(x)

    #
    #data_NCS = [[x[2],x[3]] for x in NCS]
    #data_ASU = [[x[2],x[3]] for x in ASU]
    #data_issues = [[x[2],x[3]] for x in data_files if (abs((x[1]-x[3])/x[1])>0.5)]
    ##data_issues = [[x[2],x[3]] for x in NCS if (abs((x[1]-x[3])/x[1])>0.5)]
    #plot_data4(data_NCS,data_ASU,data_issues,plt_lim=0.7)
    #
    data_NCS = [[x[2]-x[1],x[3]-x[1]] for x in NCS]
    data_ASU = [[x[2]-x[1],x[3]-x[1]] for x in ASU]
    data_issues = []
    plot_data4(data_NCS,data_ASU,data_issues,plt_lim=0.4)
    #
    #data_NCS = [[x[1],x[3]] for x in NCS]
    #data_ASU = [[x[1],x[3]] for x in ASU]
    #data_issues = []
    #plot_data4(data_NCS,data_ASU,data_issues,plt_lim=0.65)
    #
    #data = [[(x[1]-x[2])/x[2]),x[2:]] for x in data_files]
    #plot_data3(data)
    #
    #data = [[(x[1]-x[2])/x[2],x[0],x[-1]] for x in data_files]
    #data.sort()
    #print data[-20:]
    #print [x[1] for x in data[:20]]
    #
    print '*'*60
    print 'total files in plot: {}'.format(len(data_files))
    print 'number of files with >50% from published: {}'.format(len(data_issues))
    print 'Out of those, {} are from files with NCS'.format(len([x for x in NCS if (abs((x[1]-x[3])/x[1])>0.5)]))

    print '*'*60
    # collect interesting files
    file_list = [x for x in data_files if x[3]-x[2] > 0.0]
    #file_list = [x[0] for x in data_files if (x[2] == x[3])]
    print len(file_list)
    for x in file_list:
        print x
    #
    print '*'*60
    tmp = [x for x in NCS if (abs((x[1]-x[3])/x[1])>0.5)]
    print 'number of files with r-work value issues: {}'.format(len(tmp))
    for x in tmp:
        print x
    print 'done'
