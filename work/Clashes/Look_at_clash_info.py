from __future__ import division
import numpy as np
from numpy.random import rand
#import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import cPickle as pickle
import pylab as plb
import os, sys


def run():
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

    
    pdb_clash_scores.sort()
    pdb_clash_score_and_name.sort()
    print len(pdb_clash_score_and_name)
    print_list(pdb_clash_score_and_name[-100:], 5)
    print pdb_clash_score_and_name[:5]
    return pdb_clash_score_and_name

def print_list(l,n):
    '''print list l with n items in a raw''' 
    x = len(l) % n
    l.extend(['',]*x)
    for i in range(len(l)//n):
        s = i*n
        e = s + n
        print l[s:e]
    


def plot_data(pdb_clash_score_and_name):

    def onpick3(event):
        ind = event.ind
        i = ind[0]
        print '{0}  Clash Score: {1:.4f}'.format(pdb_clash_score_and_name[i][1],pdb_clash_score_and_name[i][0])
        #print 'onpick3 scatter:', ind, npy.take(x, ind), npy.take(y, ind)

    gr = 1.61803398875
    h = 10
    w = gr*h
    d = 0.03
    fig = plt.figure(figsize=(w,h))
    plt.subplots_adjust(left=d, right=1-d, top=1-d, bottom=d)
    ax1 = fig.add_subplot(111)
    #col = ax1.scatter(z, y, 100*s, c, picker=True)
    l = len(pdb_clash_score_and_name)
    x = xrange(1,l+1)
    s = 0.5 + abs(rand(l) - 1)
    c = rand(l)
    data = [d[0] for d in pdb_clash_score_and_name]
    col = ax1.scatter(x, data, 100*s, c, picker=True)
    #fig.savefig('pscoll.eps')
    fig.canvas.mpl_connect('pick_event',onpick3)
    fig.set_size_inches(w,h)
    ax1.set_ylim([0,max(data)+100])
    ax1.set_xlim([0,l+1000])
    plt.show()

if __name__=='__main__':
    pdb_clash_score_and_name = run()
    plot_data(pdb_clash_score_and_name)
    print 'done'
