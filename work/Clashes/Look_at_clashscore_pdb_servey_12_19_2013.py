from __future__ import division
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
from numpy.random import rand
import matplotlib.pyplot as plt
import cPickle as pickle
import pylab as plb
import os, sys


def get_data():
    ''' () -> dict,dict
    Read results of clash score servey for PROBE clashscore and restraints manager nonbonded clashscore
    
    c:\Phenix\Dev\Work\work\Clashes\Data\clashscore_data_dict
    
    clashscore_data_dict[pdb_file] = [clashscore_all_clashes,
                                    clashscore_simple,
                                    clashscore_only_sym_op,
                                    clashscore_solvent_solvent
                                    probe_clashscore_Ovdw_14,
                                    probe_clashscore_Ovdw_152,
                                    pdb_year]
    
    experiment_type_dict: experiment_type_dict[experiment_type] = list of pdb_file_name
    
    Returns:
    clashscore_data_dict, experiment_type_dict
    '''
    datapath = os.path.realpath('c:\Phenix\Dev\Work\work\Clashes\Data')
    clashscore_data_dict = pickle.load(open(os.path.join(datapath,'clashscore_data_dict'),'r'))
  
    experiment_dict_file = 'experiment_type_to_files_dict'		# source for experiment_type_dict
    experiment_type_dict = pickle.load(open(os.path.join(datapath,experiment_dict_file),'r'))
    
    return experiment_type_dict,clashscore_data_dict


def plot_data(experiment_type_dict,clashscore_data_dict,year_limit=1950):
    '''
    data_dict[pdb_file_name] = [total_nb_clashscore,without_sym_nb_clashscore,clashscore_probe]
    experiment_type_dict[experiment_type] = list of pdb_file_name
    pdb_file_list is a list of all the keys of data_dict
    
    ploting for each pdb file the three clashscores
    1) clashscore_probe (y1)
    2) without_sym_nb_clashscore (y2)
    3) total_nb_clashscore (y3)
    
    before printing, sort by clashscore_probe
    
    clashscore_data_dict records numbers:
    0: clashscore_all_clashes
    1: clashscore_simple,
    2: clashscore_only_sym_op,
    3: clashscore_solvent_solvent
    4: probe_clashscore_Ovdw_14,
    5: probe_clashscore_Ovdw_152,
    6: pdb_year
    '''
    
    for k in experiment_type_dict:
        # create a list with the same color for all points with the same experiment type
        pdb_files = experiment_type_dict[k]
        # collect data
        # data = [[[total_nb_clashscore,without_sym_nb_clashscore,clashscore_probe],file_name],[...],...]
        data = [[clashscore_data_dict[key],key] for key in pdb_files if clashscore_data_dict.has_key(key) and \
                (clashscore_data_dict[key][6]>year_limit)]
        # sort data by clashscore_probe 1.52
        #data.sort(key=lambda x:x[0][5])
        # sort data by clashscore_probe 1.4
        data.sort(key=lambda x:x[0][4])
        
        x = range(1,len(data)+1)
        # clashscore_probe
        #y1 = [d[0][5] for d in data]  # 1.52
        y1 = [d[0][4] for d in data]  # 1.4
        # nb_clashscore simle, without symetry or solvent-solvent clashes
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
    h = 12				# figure hight
    w = gr*h			# figure width
    d = 0.07			# distance between plot regon and figure edge
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
    fig.savefig('{}_score_compare_12_19_2013.png'.format(k))
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
    
def plot_scatter_plots(experiment_type_dict,clashscore_data_dict,year_limit=1950):
    '''
    scatter plots
    
    clashscore_data_dict records numbers:
    0: clashscore_all_clashes
    1: clashscore_simple,
    2: clashscore_only_sym_op,
    3: clashscore_solvent_solvent
    4: probe_clashscore_Ovdw_14,
    5: probe_clashscore_Ovdw_152,
    6: pdb_year
    '''
    for k in experiment_type_dict:
        # create a list with the same color for all points with the same experiment type
        pdb_files = experiment_type_dict[k]
        # collect data
        # data = [[[total_nb_clashscore,without_sym_nb_clashscore,clashscore_probe],file_name],[...],...]
        data = [[clashscore_data_dict[key],key] for key in pdb_files if clashscore_data_dict.has_key(key) and \
                (clashscore_data_dict[key]>year_limit)]
        # sort data by clashscore_probe
        #data.sort(key=lambda x:x[0][5]) # (O vdw 1.52)
        data.sort(key=lambda x:x[0][4]) # (O vdw 1.4)
 
        # clashscore_probe 
        y = [d[0][5] for d in data]  # (O vdw 1.52)
        y = [d[0][4] for d in data]  # (O vdw 1.4)
        # nb_clashscore simple
        x1 = [d[0][1] for d in data]
        # total_nb_clashscore
        x2 = [d[0][0] for d in data]
        # file names
        files = [d[1] for d in data]
        
        # set figure look
        #gr = 1.61803398875
        gr = 1
        h = 12	# figure hight
        w = gr*h	# figure width
        d = 0.07	# distance between plot regon and figure edge
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
        ax1.legend([sp1,sp2],['probe - nb_clashscore simple','probe - total_nb_clashscore'])
        fig.savefig('{}_scatter_plot_12_19_2013.png'.format(k))
        plt.show()     
        

    
def set_working_path():
    osType = sys.platform
    if osType.startswith('win'):
        directory_path = r'c:\Phenix\Dev\Work\work\Clashes'
    else:
        directory_path = '/net/cci-filer2/raid1/home/youval/Work/work/Clashes'
    # convert the path to python format
    directory_path = os.path.realpath(directory_path)
    os.chdir(directory_path)
    
def run():
    # locate the directory containing the log files
    set_working_path()
    experiment_type_dict,clashscore_data_dict = get_data()
    #data_dict,experiment_type_dict,pdb_file_list = get_data()
    
    plot_data(experiment_type_dict,clashscore_data_dict,year_limit=2000)
    plot_scatter_plots(experiment_type_dict,clashscore_data_dict,year_limit=2000)
    print 'Done!!!'

if __name__=='__main__':
    run()
