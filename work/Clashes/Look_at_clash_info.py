from __future__ import division
import numpy as np
from numpy.random import rand
#import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import cPickle as pickle
import pylab as plb
import os, sys


def run():
    '''
    Read results of clash score servey

    data = list([score_with_hydrogen,score_without_hydrogen]...)
    data_files = list([score_with_hydrogen,score_without_hydrogen,experment_type,file_name]...)
    data_dict[file_name] = [score_with_hydrogen,score_without_hydrogen,experment_type]
    '''
    # locate the directory containing the log files
    osType = sys.platform
    if osType.startswith('win'):
        directory_path = 'c:\Phenix\Dev\Work\work\Clashes'
    else:
        directory_path = '/net/cci-filer2/raid1/home/youval/Work/work/Clashes'
    # convert the path to python format
    directory_path = os.path.realpath(directory_path)
    os.chdir(directory_path)

    pdb_clash_scores = pickle.load(open('pdb_clash_scores','r'))
    pdb_clash_score_and_name = pickle.load(open('pdb_clash_score_and_name','r'))
    pdb_clash_score_dict = pickle.load(open('pdb_clash_score_dict','r'))
    # in the original run ELECTRON MICROSCOPE was not an option - fix that
    for i,x in enumerate(pdb_clash_score_and_name):
        if x[2]=='':
            pdb_clash_score_and_name[i][2] = 'ELECTRON MICROSCOPE'

    pdb_clash_scores.sort()
    pdb_clash_score_and_name.sort()
    print 'Total number of clash score records is: {}'.format(len(pdb_clash_score_and_name))
    print '*'*60
    #print_list(pdb_clash_score_and_name[-6:], 2)
    #print_list(pdb_clash_score_and_name[:50], 5)
    return pdb_clash_score_and_name,pdb_clash_score_dict

def print_list(l,n):
    '''print list l with n items in a raw'''
    x = len(l) % n
    l.extend(['',]*x)
    for i in range(len(l)//n):
        s = i*n
        e = s + n
        print l[s:e]



def plot_data(pdb_clash_score_and_name,by_type_dict):

    def onpick3(event):
        ind = event.ind
        i = ind[0]
        print '{0}  Clash Score: {1:.4f}'.format(pdb_clash_score_and_name[i][1],pdb_clash_score_and_name[i][0])
        #print 'onpick3 scatter:', ind, npy.take(x, ind), npy.take(y, ind)

    # set figure look
    gr = 1.61803398875
    h = 10
    w = gr*h
    d = 0.03
    fig = plt.figure(figsize=(w,h))
    plt.subplots_adjust(left=d, right=1-d, top=1-d, bottom=d)
    ax1 = fig.add_subplot(111)
    #
    #col = ax1.scatter(z, y, 100*s, c, picker=True)
    #l = len(pdb_clash_score_and_name)
    #x = xrange(1,l+1)
    # size of data point reflects the diference between the with - without hydrogen scores
    #s = [0.5+abs(x[0]-x[1]) for x in pdb_clash_score_and_name]
    #s = 0.5 + abs(rand(l) - 1)
    # set the colors. there are 6 type of experiments
    ct = [x/6 for x in range(1,7)]
    c = []
    data = []
    experiment_type = []
    for i,k in enumerate(by_type_dict):
        # create a list with the same color for all points with the same experiment type
        c.extend([ct[i]]*len(by_type_dict[k]))
        # collect data
        experiment_type.append(k)
        data.extend(by_type_dict[k])
    #data = [d[0] for d in data]
    #col = ax1.scatter(x, data, 100*s, c, picker=True)
    # build data with size and color
    data = [[i,d[0],100*(0.5+abs(d[0]-d[1])),c[i]] for i,d in enumerate(data)]
    x = [d[0] for d in data]
    y = [d[1] for d in data]
    s = [d[2] for d in data]
    c = [d[3] for d in data]
    col = ax1.scatter(x,y,s,c, picker=True)
    #fig.savefig('pscoll.eps')
    fig.canvas.mpl_connect('pick_event',onpick3)
    fig.set_size_inches(w,h)
    ax1.set_ylim([0,max(y)+100])
    ax1.set_xlim([0,x[-1]+1000])
    plt.show()

def zero_scores(clash_data):
    '''
    In ELECTRON MICROSCOPE files there is no EXPERIMENT TYPE and clash score are zero
    Remove those records
    '''
    zero_score = [x for x in clash_data if x[0:2] == [0,0]]
    zero_scores_dict = initial_dict()
    for x in zero_score:
        keys = x[2].split(',')
        for k in keys:
            zero_scores_dict[k].append(x[3])

    print 'Number of records with 0.0 clash scores'
    print '='*60
    for x in zero_scores_dict:
        print 'from {0} is: {1}'.format(x,len(zero_scores_dict[x]))
    print '*'*60

def create_by_type_dict(pdb_clash_score_and_name):
    '''(list) -> dicttionary
    sort clash score by experiment type
    '''
    by_type_dict = initial_dict()
    for x in pdb_clash_score_and_name:
        keys = x[2].split(',')
        for k in keys:
            by_type_dict[k].append([x[0],x[1],x[3]])

    print 'Experimental type breakdown'
    print '='*60
    for x in by_type_dict:
        print 'Number of from {0} is: {1}'.format(x,len(by_type_dict[x]))
    print '*'*60
    return by_type_dict

def initial_dict():
    init_dict = dict([('X-RAY DIFFRACTION',[]),
                      ('NMR',[]),
                      ('NEUTRON DIFFRACTION',[]),
                      ('ELECTRON MICROSCOPE',[]),
                      ('Other',[]),
                      ('SMALL ANGLE X-RAY SCATTERING',[])])
    return init_dict

if __name__=='__main__':
    pdb_clash_score_and_name,pdb_clash_score_dict = run()
    # Look at records with 0.0 scores
    zero_scores(pdb_clash_score_and_name)
    by_type_dict = create_by_type_dict(pdb_clash_score_and_name)
    plot_data(pdb_clash_score_and_name,by_type_dict)
    print 'done'
