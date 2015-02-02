from __future__ import division
from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
import cPickle as pickle
import numpy as np
import sys
import os


class Plot_data(object):

  def __init__(self):
    """
    When using multiple lines of similar data, collect the plot info using
    this object

    count (y value): count of files for each year
    percent (y value): 100* count/total_structure
    years (x value): year when PDB structure was published
    n_ignored (int): the number of structures not included
    nbo (int): file included are those with NBO >= nbo
    """
    self.count = None
    self.percent = None
    self.years = None
    self.n_ignored = None
    self.nbo = None

class Collecting_overlaps_data(object):

  def __init__(self,start_year=1980,end_year=2014):
    """
    if "test_clean_data" and "test_clean_dict" exist, plotting is done on
    existing data

    if "test_clean_data" and "test_clean_dict" do not exist, collect info from the
    folder containing the tests, from the queue and create them.

    the file "test_data.txt"  is a comma separated file, containing all data.
    the file "test_clean_data" is a pickled list containing only files that
    could be processed (without error of any kind).
    the file "test_clean_dict" is a pickled dictionary containing all data.

    "test_data.txt" is a comma separated file with the header:
      PDB ID,Macro molecule overlaps,Symmetry overlaps,All overlaps,
      Macro molecule overlaps per 1000 atoms,Symmetry overlaps per 1000 atoms,
      All overlaps per 1000 atoms,year model deposited in PDB,experiment type

    Args:
      start_year, end_year (int): start and end years to show on plots
    """
    self.data_file_name = 'test_data.txt'
    self.clean_dict_file_name = 'test_clean_dict'
    self.clean_data_file_name = 'test_clean_data'
    #
    self.data = []
    self.clean_data = []
    self.clean_data_dict = {}
    self.queue_data_path = "/net/cci/youval/work/work/Clashes/queue_clash_compare"
    # filtered data
    self.plot_data_list = []
    # total data
    self.n_total = None
    self.n_total_dict = {}
    self.years = None
    self.start_year = start_year
    self.end_year = end_year
    # collecting years for structures with x[i] > min_nbo
    self.nbo_per_1000_atoms = []
    self.sym = True

  def set_working_path(self):
    """ Set working and data folders """
    p = r'C:\Users\Youval\Google Drive\Documents\LBNL\phenix\publications'
    p += r'\news letter\clash score\related code'
    p += r'\pdb_overlap_scan'
    osType = sys.platform
    if osType.startswith('win'):
      assert os.path.isdir(p)
      self.working_path = p
      print 'data path: ',p
    else:
      path = '/net/cci/youval/work/work/Clashes'
      assert os.path.isdir(path)
      self.working_path = path
      print 'data path: ',path
    os.chdir(self.working_path)

  def get_test_data(self):
    """ collect data from pdb scan or existing data files  """
    have_data = os.path.isfile(self.data_file_name)
    have_data &= os.path.isfile(self.clean_dict_file_name)
    have_data &= os.path.isfile(self.clean_data_file_name)
    if have_data:
      print 'using existing data files'
      self.data = self.read_csv_data()
      self.clean_data = pickle.load(open(self.clean_data_file_name,'r'))
      self.clean_data = [x for x in self.clean_data if x[1] >= 0 ]
      self.data_dict = pickle.load(open(self.clean_dict_file_name,'r'))
      print "Number of good files: ",len(self.clean_data)
      print "Total number files processed: ",len(self.data)
    else:
      print 'getting new data from {}'.format(self.queue_data_path)
      # check if data folder exist
      if os.path.isdir(self.queue_data_path):
        # Read files in directory_path
        files = os.listdir(self.queue_data_path)
        # collect only the files that starts with log_
        files = [x for x in files if x.startswith('log_')]
        print "Number of log files: ",len(files)
        for fn in files:
          d = open(os.path.join(self.queue_data_path, fn), "r").readlines()
          if d:
            data = format_data_types(d[0])
          else:
            data = []
          if not ((len(d)==1) and (len(data) == 9)):
            # Some issue with results
            pdb_id = fn[-4:]
            data  = [pdb_id] + ([-9] * 8)
          self.data.append(data)
        # clean data, collect good data
        print 'Total number data records: {}'.format(len(self.data))
        f = open(self.data_file_name,'w')
        out_str = ['{}'] * 9
        out_str = ','.join(out_str) + '\n'
        for d in self.data:
          pdb_id = d[0]
          f.write(out_str.format(*d))
          if d[1] >= 0:
            self.clean_data.append(d)
            self.clean_data_dict[pdb_id] = d
        f.close()
        print "Number of good records: ",len(self.clean_data)
        pickle.dump(self.clean_data, open(self.clean_data_file_name,'w'))
        pickle.dump(self.clean_data_dict, open(self.clean_dict_file_name,'w'))

  def read_csv_data(self):
    """ read the data from csv text file """
    data = open(self.data_file_name,'r').read().splitlines()
    data_list = []
    for l in data:
      if (not l.startswith('#')) and (not l.startswith('PDB ID')):
        d = format_data_types(l)
        data_list.append(d)
    return data_list

  def prepare_data_for_plotting(self,sym=True):
    """
    Process data for plotting.
    nonbonded overlaps (NBO) > min_nbo vs.
    year pdb structure was submitted, starting at 1980 till 2014

    Build self.plot_data_list

    Args:
      sym (bool): when True plot symmetry NBO when False use all NBO
    """
    self.sym = sym
    plt.close('all')
    if not self.clean_data:
      raise IOError('No Data to plot')
    # get data year: x[7], sym NBO: x[2], all NBO: x[3]
    # get data year: x[7], sym NBO per 1000 atoms: x[5], all NBO: x[3]
    assert len(self.clean_data[0]) == 9
    if sym:
      # i = 2 # NBO
      i = 5 # NBO per 1000 atoms
    else:
      # i = 3 # NBO
      i = 6 # NBO per 1000 atoms
    # collecting years for structures with x[i] > min_nbo
    # nbo_list = [10,20,50,100]
    nbo_per_1000_atoms = [0,3,6,9,15]
    data = [x[7] for x in self.clean_data
            if (x[7] >= self.start_year) and (x[7] <= self.end_year)]
    n_bins = np.arange(min(data)-0.5,max(data)+1.5,1)
    years =  range(min(data),max(data)+1,1)
    # get the total of structure for each year
    n_total, all_years, _ = plt.hist(data, bins=n_bins)
    n_total_dict = {y:n for y,n in zip(years,n_total)}
    plt.clf()
    n_all = len(self.clean_data)
    # for min_nbo in nbo_list:
    for min_nbo in nbo_per_1000_atoms:
      pd = Plot_data()
      pd.nbo = min_nbo
      data = [x[7] for x in self.clean_data if x[i] > min_nbo]
      pd.n_ignored = len(self.clean_data)-len(data)
      # use histogram to get the data point for each nbo
      pd.count, _, _ = plt.hist(data, bins=n_bins)
      plt.clf()
      pd.years = np.array(years)
      assert len(pd.count) == len(years)
      # collect the total number of structure deposited for present years
      tmp = np.array([n_total_dict[y] for y in pd.years])
      tmp = np.array([(x,y) for x,y in zip(pd.count,tmp) if y>0])
      pd.percent = 100*tmp[:,0]/tmp[:,1]
      pd.count = tmp[:,0]
      assert len(pd.years) == len(pd.percent)
      # filtered data
      self.plot_data_list.append(pd)
      # total data
      self.n_total = n_total
      self.years = years
      self.nbo_per_1000_atoms = nbo_per_1000_atoms
      self.n_total_dict = n_total_dict

  def nbo_vs_year(self):
    """
    Plot percent of structures, with different NBO per 1000 atoms levels,
    from "good" pdb structures (all PDB files with a single model, no unknown
    atom types and good CRYST1 records) VS. year

    Second sub plot: the total of "good" structures deposited VS. year
    """
    plt.close('all')
    # figure parameters
    # plt.ion() # enables interactive mode
    max_y = 105
    fontsize = 20
    fig = plt.figure(figsize=(8,10))
    gs = GridSpec(2,1,height_ratios=[2,1])
    # first subplot
    ax1 = plt.subplot(gs[0,0])
    ax2 = plt.subplot(gs[1,0])
    lines = []
    line_type = ['.:','.-','.--']
    n = len(line_type)
    for i,pd in enumerate(self.plot_data_list):
      lt = line_type[i%n]
      l, = ax1.plot(pd.years,pd.percent,lt)
      lines.append(l)
    ax1.set_ylabel('Percent of PDB structures',fontsize=fontsize)
    ax1.text(min(self.years)+0.5,max_y-4,'a.',fontsize=fontsize)
    ax1.tick_params(axis='both',labelsize=fontsize - 2)
    ax1.axes.get_xaxis().set_visible(False)
    ax1.set_yticks([5,10,40,70,100])
    ax1.set_ylim([0,max_y])
    ax1.set_xlim([self.start_year,self.end_year])
    # legend
    labels = ['NBO per 1000 atom > {}']*len(self.nbo_per_1000_atoms)
    labels = [x.format(y) for x,y in zip(labels,self.nbo_per_1000_atoms)]
    if self.sym:
      legend_pos = [0.96,0.70]
    else:
      legend_pos = [0.54,0.30]
    ax1.legend(
      lines,labels,
      bbox_to_anchor=legend_pos,
      loc=1,borderaxespad=0.0)
    # Second subplot
    ax2.plot(self.years,self.n_total,'.:g')
    ax2.set_xlim([self.start_year,self.end_year])
    ax2.set_xlabel('Year',fontsize=fontsize)
    ax2.set_ylabel('Number of structures',fontsize=fontsize)
    ax2.text(min(self.years)+0.5,max(self.n_total)-5,'b.',fontsize=fontsize)
    ax2.tick_params(axis='both',labelsize=fontsize - 2)
    ax2.set_xticks([self.start_year,1990,2000,self.end_year])
    ax2.set_yscale('log')
    ax2.set_yticks([10,100,1000])
    #
    gs.tight_layout(fig)
    gs.update(hspace=0)
    s = 'all'*(not self.sym) + 'sym'*self.sym
    fig_name = 'nbo_vs_year_{}.png'.format(s)
    plt.savefig(fig_name)
    fig.show()

  def show_values(self,year=2014):
    """ Show the plot numbers for a particular year """
    print '--------------------------'
    print 'Showing info for {}'.format(year)
    print '--------------------------'
    msg = 'total number of structure PDB files with a single model, \n'
    msg += 'no unknown atom types and good CRYST1 records: {}\n'
    print msg.format(int(self.n_total_dict[year]))
    # find and collect number of structure. for each NBO limit,
    # for the desired year
    data = self.plot_data_list[0]
    i = list(data.years).index(year)
    n = [x.count[i] for x in self.plot_data_list]
    p = [x.percent[i] for x in self.plot_data_list]
    #
    s = ['NBO per 1000 atom > {:3}, , # of structure: {:5}, percent: {:5.1f}']
    s = s*len(self.nbo_per_1000_atoms)
    s = [x.format(y,int(z),k) for x,y,z,k in zip(s,self.nbo_per_1000_atoms,n,p)]
    for l in s:
      print l
    print

def format_data_types(s):
  """
  apply the correct data type to each value in the list created from a comma
  separated sting "s"

  x1: PDB ID (string)
  x2: Macro molecule overlaps (int)
  x3: Symmetry overlaps (int)
  x4: All overlaps (int)
  x5: Macro molecule overlaps per 1000 atoms (float)
  x6: Symmetry overlaps per 1000 atoms (float)
  x7: All overlaps per 1000 atoms (float)
  x8: year model deposited in PDB (int)
  x9: experiment type (string)
  """
  d = [x.strip() for x in s.split(',')]
  if len(d) == 9:
    # integer values
    for i in [1,2,3,7]:
      d[i] = int(d[i])
    # float values
    for i in [4,5,6]:
      d[i] = round(float(d[i]),1)
    return d
  else:
    return None

def run():
  print 'Start'
  test_results = Collecting_overlaps_data()
  test_results.set_working_path()
  test_results.get_test_data()
  test_results.prepare_data_for_plotting(sym=False)
  test_results.nbo_vs_year()
  test_results.show_values(year=2014)


  print 'Done...'


if __name__=='__main__':
  run()
