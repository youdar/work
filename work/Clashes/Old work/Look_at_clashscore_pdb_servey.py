from __future__ import division
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
from numpy.random import rand
import matplotlib.pyplot as plt
import cPickle as pickle
import pylab as plb
import os, sys


def get_data():
    '''
    Read results of clash score servey for PROBE clashscore and restraints manager nonbonded clashscore
    
    
    c:\Phenix\Dev\Work\work\Clashes\Data\clashscore_compare_ready_set-12-5-13_dict
    data_dict[pdb_file_name] = [total_nb_clashscore,without_sym_nb_clashscore,clashscore_probe]
    experiment_type_dict[experiment_type] = list of pdb_file_name
    
    >>>experiment_type_dict['NMR'][:10]
    ['103d', '124d', '141d', '142d', '169d', '175d', '1a1d', '1ac3', '1al5', '1anp']
    >>>data_dict['142d']
    [0.0, 0.0, 0.0]

    
    pdb_clash_scores = list([score_with_hydrogen,score_without_hydrogen]...)
    pdb_clash_score_and_name = list([score_with_hydrogen,score_without_hydrogen,experment_type,file_name]...)
    pdb_clash_score_dict[file_name] = [score_with_hydrogen,score_without_hydrogen,experment_type]
    '''
    # Get data
    path = os.path.realpath('c:\Phenix\Dev\Work\work\Clashes\Data')
    #f = os.path.join(path,'clashscore_compare_reduce_12_6_2013_dict')  # Probe O vdw is 1.4
    f = os.path.join(path,'clashscore_compare_reduce_12_11_2013_dict') # Probe O vdw is 1.52
    data_dict = pickle.load(open(f,'r'))
    f = os.path.join(path,'experiment_type_to_files_dict')
    experiment_type_dict = pickle.load(open(f,'r'))
    # Collect all files that we compared 
    pdb_file_list = [key for key in data_dict]
    return data_dict,experiment_type_dict,pdb_file_list

def plot_data(data_dict,experiment_type_dict,pdb_file_list):
    '''
    data_dict[pdb_file_name] = [total_nb_clashscore,without_sym_nb_clashscore,clashscore_probe]
    experiment_type_dict[experiment_type] = list of pdb_file_name
    pdb_file_list is a list of all the keys of data_dict
    
    ploting for each pdb file the three clashscores
    1) clashscore_probe (y1)
    2) without_sym_nb_clashscore (y2)
    3) total_nb_clashscore (y3)
    
    before printing, sort by clashscore_probe
    '''
    
    for k in experiment_type_dict:
        # create a list with the same color for all points with the same experiment type
        pdb_files = experiment_type_dict[k]
        # collect data
        # data = [[[total_nb_clashscore,without_sym_nb_clashscore,clashscore_probe],file_name],[...],...]
        data = [[data_dict[key],key] for key in pdb_files if data_dict.has_key(key)]
        # sort data by clashscore_probe
        data.sort(key=lambda x:x[0][2])
        
        x = range(1,len(data)+1)
        # clashscore_probe
        y1 = [d[0][2] for d in data]
        # without_sym_nb_clashscore
        y2 = [d[0][1] for d in data]
        # total_nb_clashscore
        y3 = [d[0][0] for d in data]
        # file names
        files = [d[1] for d in data]

        plot_experiment(x,y1,y2,y3,k,files)  
        #hist_both_clash_scores(y,y2,k)
    

def plot_experiment(x,y1,y2,y3,k,files):
    '''
    plot_experiment(x,y,s,c,k,data):
    plot a sub plot for an experiment type
    x: enumerating data points
    y1: clashscore_probe
    y2: without_sym_nb_clashscore
    y3: total_nb_clashscore
    k: pdb file experiment type
    files: pdb file names
    '''    
    def line_picker(p1, mouseevent):
        """
        find the points within a certain distance from the mouseclick in
        data coords and attach some extra attributes, pickx and picky
        which are the data points that were picked
        """
        if mouseevent.xdata is None: return False, dict()
        xdata = p1.get_xdata()
        ydata = p1.get_ydata()
        maxd = 0.5
        d = np.sqrt((xdata-mouseevent.xdata)**2.)

        ind = np.nonzero(np.less_equal(d, maxd))
        if len(ind):
            pickx = np.take(xdata, ind)
            picky = np.take(ydata, ind)
            props = dict(ind=ind, pickx=pickx, picky=picky)
            i =  pickx-1
            if i.size!=0:
                print '*'*50
                try:
                    print 'PDB file {0}\nExperiment type: {1}'.format(files[i],k)
                except TypeError:
                    pass
                print 'clashscore_probe: {0:.4f}\nwithout_sym_nb_clashscore: {1:.4f}\ntotal_nb_clashscore: {2:.4f}'.format(y1[i],y2[i],y3[i])
                return True, {}
        else:
            return False, dict()

    # set figure look
    #gr = 1.61803398875 	# Golden ratio
    gr = 1
    h = 10				# figure hight
    w = gr*h			# figure width
    d = 0.05			# distance between plot regon and figure edge
    fig = plt.figure(figsize=(w,h))
    plt.subplots_adjust(left=d, right=1-d, top=1-d, bottom=d)
    ax1 = fig.add_subplot(111)
    # set scattering plot and allow interactinve selection of points on plot
    p2, = ax1.plot(x,y2,'.y')
    p3, = ax1.plot(x,y3,'.b')
    p1, = ax1.plot(x,y1,'.r', picker=line_picker)
    fig.set_size_inches(w,h)
    ax1.legend([p1,p2,p3],['clashscore_probe','without_sym_nb_clashscore','total_nb_clashscore'])
    #
    maxy = 100
    maxs = 0.5
    ax1.set_ylim([-maxy*.01,maxy+maxs])
    ax1.set_xlim([-x[-1]*0.01,x[-1]*1.01+maxs])
    #
    plt.title(k)
    plt.ylabel('clashscore and nb_clashscore')
    #fig.savefig('{}.eps'.format(k))
    fig.savefig('{}_score_compare.png'.format(k))
    plt.show()     


def hist_both_clash_scores(x,y,k):
    '''
    x: clash score without pdb hydrogens as y (keep_hydrogens=False)
    y: clash score with pdb hydrogens as y
    k: Experiment type
    '''
    
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
    
def plot_scatter_plots(data_dict,experiment_type_dict,pdb_file_list):
    '''
    scatter plots
    '''
    for k in experiment_type_dict:
        # create a list with the same color for all points with the same experiment type
        pdb_files = experiment_type_dict[k]
        # collect data
        # data = [[[total_nb_clashscore,without_sym_nb_clashscore,clashscore_probe],file_name],[...],...]
        data = [[data_dict[key],key] for key in pdb_files if data_dict.has_key(key)]
        # sort data by clashscore_probe
        data.sort(key=lambda x:x[0][2])
 
        # clashscore_probe
        y = [d[0][2] for d in data]
        # without_sym_nb_clashscore
        x1 = [d[0][1] for d in data]
        # total_nb_clashscore
        x2 = [d[0][0] for d in data]
        # file names
        files = [d[1] for d in data]
        
        # set figure look
        #gr = 1.61803398875
        gr = 1
        h = 10	# figure hight
        w = gr*h	# figure width
        d = 0.05	# distance between plot regon and figure edge
        fig = plt.figure(figsize=(w,h))
        plt.subplots_adjust(left=d, right=1-d, top=1-d, bottom=d)
        ax1 = fig.add_subplot(111)
        # set scattering plot and allow interactinve selection of points on plot
        sp1 = ax1.scatter(x1,y,c='y')
        sp2 = ax1.scatter(x2,y)
        fig.set_size_inches(w,h)
        #
        maxy = 100
        maxs = 0.5
        ax1.set_ylim([0,100])
        ax1.set_xlim([0,100])
        #
        plt.title(k)
        plt.ylabel('clashscore')
        plt.xlabel('nb_clashscore')
        fig.savefig('{}_scatter_plot.png'.format(k))
        plt.show()     
        

    
def set_working_path():
    osType = sys.platform
    if osType.startswith('win'):
        directory_path = 'c:\Phenix\Dev\Work\work\Clashes'
    else:
        directory_path = '/net/cci-filer2/raid1/home/youval/Work/work/Clashes'
    # convert the path to python format
    directory_path = os.path.realpath(directory_path)
    os.chdir(directory_path)
    
def run():
    # locate the directory containing the log files
    set_working_path()
    data_dict,experiment_type_dict,pdb_file_list = get_data()
    #plot_data(data_dict,experiment_type_dict,pdb_file_list)
    plot_scatter_plots(data_dict,experiment_type_dict,pdb_file_list)
    print 'Done!!!'

if __name__=='__main__':
    run()
