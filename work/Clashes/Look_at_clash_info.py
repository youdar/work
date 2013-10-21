from __future__ import division
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
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
    
    for k in by_type_dict:
        # create a list with the same color for all points with the same experiment type
        data = by_type_dict[k]
        
        #c = np.ones(len(data))*0.647933889333
        # build data with size and color
        #datalist = [[i,d[0],(d[0]-d[1])] for i,d in enumerate(data)]
        x = range(1,len(data)+1)
        #x = [d[0] for d in datalist]
        y = [d[1] for d in data]	# use clash score without pdb hydrogens as y
        y2 = [d[0] for d in data]	# use clash score with pdb hydrogens as y
        # make the size of the points on the plot relative to the difference in the clash scores
        s = [50 + 5*abs(d[1]-d[0]) for d in data]
        # The color of points where both clash scores are the same
        c = ['y',]*len(data)
        # Color the data points in a different colors
        for i in range(len(data)):
            if data[i][0]>data[i][1]: c[i] = 'b'
            elif data[i][0]<data[i][1]: c[i] = 'r'
        #c = rand(len(data))
        plot_experiment(x,y,s,c,k,data)  
        #hist_both_clash_scores(y,y2,k)
    

def plot_experiment(x,y,s,c,k,data):
    '''
    plot a sub plot for an experiment type
    x: enumerating data points
    y: clash score with hydrogen
    s: size the data point, related to the difference between with/without hydrogen clash scores
    c: data point color
    k: pdb file experiment type
    '''
    def onpick3(event):
        ind = event.ind
        i = ind[0]
        print '*'*50
        print 'PDB file {0}     Experiment type: {1}'.format(data[i][2],k)
        print 'Clash score with hydrogen kept: {0:.4f}    without hydrogen: {1:.4f}'.format(data[i][0],data[i][1])
        print c[i]

    # set figure look
    gr = 1.61803398875
    h = 10	# figure hight
    w = gr*h	# figure width
    d = 0.05	# distance between plot regon and figure edge
    fig = plt.figure(figsize=(w,h))
    plt.subplots_adjust(left=d, right=1-d, top=1-d, bottom=d)
    ax1 = fig.add_subplot(111)
    # set scattering plot and allow interactinve selection of points on plot
    col = ax1.scatter(x,y,s,c=c, picker=True)
    fig.canvas.mpl_connect('pick_event',onpick3)
    fig.set_size_inches(w,h)
    #
    maxy = max(y)
    maxs = max(s)/100
    ax1.set_ylim([-maxy*.01,maxy+maxs])
    ax1.set_xlim([-x[-1]*0.01,x[-1]*1.01+maxs])
    #
    plt.title(k)
    delta_score = [abs(i[0]-i[1]) for i in data]
    minscore = min(delta_score)
    maxscore = max(delta_score)
    text1 = 'Number of data points: {0}\nMin score difference: {1}\nMax score difference: {2}\n\n'.format(x[-1],minscore,maxscore)
    text2 = 'Blue: Score excluding H is lower\nRed: Score including PDB H is lower\nYellow: The same '
    plt.text(x[-1]*0.1,maxy*.65, text1+text2,fontsize=16)
    plt.ylabel('Clash score - Excluding hydrogens in input file')
    fig.savefig('pscoll.eps')
    plt.show()     


def hist_both_clash_scores(x,y,k):
    
    # set figure look
    
    h = 11	# figure hight
    w = 11	# figure width
    d = 0.05	# distance between plot regon and figure edge
    plt.figtext(0,0,k, fontsize=16)
    fig, axScatter = plt.subplots(figsize=(w,h))

    plt.subplots_adjust(left=d, right=1-d, top=1-d, bottom=d)

    # the scatter plot:
    axScatter.scatter(x, y)
    axScatter.set_aspect(1.)
    
    # create new axes on the right and on the top of the current axes
    # The first argument of the new_vertical(new_horizontal) method is
    # the height (width) of the axes to be created in inches.
    divider = make_axes_locatable(axScatter)
    binlim = 200
    axHistx = divider.append_axes("top", 3, pad=1, sharex=axScatter, xlabel='Clash score without PDB Hydrogen', xlim=[0,binlim])
    axHisty = divider.append_axes("right", 3, pad=1, sharey=axScatter, ylabel='Clash score with PDB Hydrogen', ylim=[0,binlim])
    plt.annotate
    
    #bins = np.arange(-lim, lim + binwidth, binwidth)
    bins = 40
    axHistx.hist(x, bins=bins)
    axHisty.hist(y, bins=bins, orientation='horizontal')
    plt.figtext(0.4,0.97,k, fontsize=16)
    plt.draw()
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

    print 'Number of records with 0.0 clash scores: {}'.format(len(zero_score))
    print '='*60
    for x in zero_scores_dict:
        print '{0:30} : {1:4}'.format(x,len(zero_scores_dict[x]))
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
        print ' {0:30} : {1:4}'.format(x,len(by_type_dict[x]))
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
