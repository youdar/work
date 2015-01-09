from __future__ import division
from scipy.stats import linregress
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt
import cPickle as pickle
import pylab as plb
import os, sys

class Collecting_clashscore_data(object):

  def __init__(self):
    """
    Collect and analyse clashscore test data

    good test result should have only one line
      x1::x2::x3::x4::x5
      Where:
        x1: pdb_file_name
        x2: macro_molecule_cctbx_clashscore
        x3: symmetry_cctbx_clashscore
        x4: total_cctbx_clashscore
        x5: probe_clashscore

    if test results contains only -1: PDB file was not found
    if test results contains only -2: there where issues getting the
      clashscore and proper structure factor file was not available
    if test results start with "Sorry" set values to -3
    """
    self.current_dir = os.getcwd()
    self.working_path = ''
    self.data_file_name = 'test_data'
    self.clean_data_file_name = 'test_clean_data'
    self.data_dict_file_name = 'test_data_dict'
    self.files_with_problem_file_name = 'files_with_issues.txt'
    self.data = []
    self.clean_data = []
    self.pdb_id_with_issues = []
    self.data_dict = {}
    self.queue_data_path = "/net/cci/youval/Work/work/Clashes/queue_clash_compare"

  def set_working_path(self):
    """ Set working and data folders """
    osType = sys.platform
    if osType.startswith('win'):
        self.working_path = r'c:\Phenix\Dev\Work\work\Clashes'
    else:
        path = '/net/cci-filer2/raid1/home/youval/Work/work/Clashes'
        self.working_path = path
    # convert the path to python format
    self.working_path = os.path.realpath(self.working_path)
    os.chdir(self.working_path)

  def change_to_original_path(self):
    """ change current directory to original one """
    os.chdir(self.current_dir)

  def get_test_data(self):
    """
    If test_data.txt and test_data_dict exit, use them, otherwise create them
    """
    have_data = os.path.isfile(self.data_file_name)
    have_data &= os.path.isfile(self.data_dict_file_name)
    if have_data:
      print 'using existing data files'
      self.data = pickle.load(open(self.data_file_name,'r'))
      self.clean_data = pickle.load(open(self.clean_data_file_name,'r'))
      self.data_dict = pickle.load(open(self.data_dict_file_name,'r'))
      self.pdb_id_with_issues = open(
        self.files_with_problem_file_name,'r').read().splitlines()
      print "Number of good file with issues: ",len(self.pdb_id_with_issues)
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
            raw_data = d[0].split('::')
          else:
            raw_data = []
          if (len(d)==1) and (len(raw_data) == 5):
            pdb_id = raw_data[0]
            data = [round(float(x),1) for x in raw_data[1:]]
            data = [pdb_id] + data
          else:
            # Some issue with results
            pdb_id = fn[-4:]
            data  = [pdb_id,-3,-3,-3,-3]
          self.data.append(data)
        # clean data, collect good data
        print 'Total number data records: {}'.format(len(self.data))
        for d in self.data:
          pdb_id = d[0]
          if (d[1] != -1) and (d[1] != -2) and (d[1] != -3):
            self.clean_data.append(d)
            self.data_dict[pdb_id] = d
          else:
            self.pdb_id_with_issues.append(pdb_id)
        print "Number of good records: ",len(self.clean_data)
        pickle.dump(self.data, open(self.data_file_name,'w'))
        pickle.dump(self.clean_data, open(self.clean_data_file_name,'w'))
        pickle.dump(self.data_dict, open(self.data_dict_file_name,'w'))
        pdb_id_with_issues = '\n'.join(self.pdb_id_with_issues)
        open(self.files_with_problem_file_name,'w').write(pdb_id_with_issues)
    print 'Number of good data points: {}'.format(len(self.clean_data))
    print '-'*50

  def plot_reference(self,ignore_delta=100):
    """
    Compare CCBTX simple non-bonded clashscore to PROBE clashscore

    Args:
      ignore_delta (float): ignore outlier points where
        abs(cctbx_score - probe_score) > ignore_delta
    """
    if self.clean_data:
      # figure limits
      max_x = 100
      max_y = 100
      fontsize = 20
      # get data
      cctbx_prob = [(x[5],x[1]) for x in self.clean_data
                    if abs(x[1]-x[5])<ignore_delta]
      n_ignored = len(self.clean_data)-len(cctbx_prob)
      n_all = len(self.clean_data)
      print 'Number of data points ignored: ',n_ignored
      cctbx_prob.sort()
      # cctbx
      cctbx_score = [x[1] for x in cctbx_prob]
      probe_score = [x[0] for x in cctbx_prob]
      # Get linear fitting parameters
      x_fit = [0,max_x]
      # m,b = plb.polyfit(cctbx_score, probe_score, 1)
      m,b ,r_value,_,_ = linregress(cctbx_score, probe_score)
      print 'r-value: ',r_value
      y_fit = [b,m * max_x + b]
      # gr = 1.61803398875 	# Golden ratio
      gr = 1
      h = 12				# figure height
      w = gr*h			# figure width
      fig = plt.figure(figsize=(w,h))
      plt.plot(cctbx_score,probe_score,'.b',x_fit,y_fit,'y',linewidth=2)
      plt.xticks(fontsize=fontsize)
      plt.yticks(fontsize=fontsize)
      plt.title(
        'CCBTX macro molecule vs PROBE non-bonded clashscore',
        fontsize=fontsize)
      plt.ylabel('PROBE clashscore',fontsize=fontsize)
      plt.xlabel('CCBTX simple non-bonded clashscore',fontsize=fontsize)
      ignore_str = 'Ignore {} of {} outlier points where\n'
      ignore_str +=  'abs(macro_mol cctbx_score - probe_score) >= {}'
      ignore_str = ignore_str.format(n_ignored,n_all,ignore_delta)
      # these are matplotlib.patch.Patch properties
      props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
      plt.text(2,max_y - 2,ignore_str,
               fontsize=fontsize,verticalalignment='top',bbox=props)
      plt.xlim([0,max_x])
      plt.ylim([0,max_y])
      #fig.savefig('{}.eps'.format(k))
      fig.savefig('cctbx_simple_vs_probe.png')
      plt.show()
      print 'Linear fit info: PROBE = {} * CCBTX + {}'.format(m,b)
      print ignore_str
      print '-'*50
    else:
      print 'Load data before attempting to plot'

  def plot_total(self,ignore_delta=100):
    """
    Compare CCBTX simple non-bonded + symmetry related clashscore
    to PROBE clashscore

    Args:
      ignore_delta (float): ignore outlier points where
        abs(cctbx_score - probe_score) > ignore_delta
    """
    if self.clean_data:
      # figure limits
      max_x = 100
      max_y = 100
      fontsize = 20
      # get data
      assert len(self.clean_data[0]) == 6
      cctbx_prob = [(x[5],x[1]+x[2]) for x in self.clean_data
                    if abs(x[1]+x[2]-x[5])<ignore_delta]
      n_ignored = len(self.clean_data)-len(cctbx_prob)
      n_all = len(self.clean_data)
      print 'Number of data points ignored: ',n_ignored
      cctbx_prob.sort()
      # cctbx
      cctbx_score = [x[1] for x in cctbx_prob]
      probe_score = [x[0] for x in cctbx_prob]
      # Get linear fitting parameters
      x_fit = [0,max_x]
      # m,b = plb.polyfit(cctbx_score, probe_score, 1)
      m,b ,r_value,_,_ = linregress(cctbx_score, probe_score)
      print 'r-value: ',r_value
      y_fit = [b,m * max_x + b]
      #gr = 1.61803398875 	# Golden ratio
      gr = 1
      h = 12				# figure height
      w = gr*h			# figure width
      fig = plt.figure(figsize=(w,h))
      plt.plot(cctbx_score,probe_score,'.b',x_fit,y_fit,'y',linewidth=2)
      plt.xticks(fontsize=fontsize)
      plt.yticks(fontsize=fontsize)
      plt.title(
        'CCBTX simple non-bonded + symmetry related vs PROBE',fontsize=fontsize)
      plt.ylabel('PROBE clashscore',fontsize=fontsize)
      plt.xlabel(
        'CCBTX simple non-bonded + symmetry related clashscore',
        fontsize=fontsize)
      ignore_str = 'Ignore {} of {} outlier points where\n'
      ignore_str +=  'abs(simple cctbx_score - probe_score) >= {}'
      ignore_str = ignore_str.format(n_ignored,n_all,ignore_delta)
      # these are matplotlib.patch.Patch properties
      props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
      plt.text(2,max_y - 2,ignore_str,
               fontsize=fontsize,verticalalignment='top',bbox=props)
      plt.xlim([0,max_x])
      plt.ylim([0,max_y])
      #fig.savefig('{}.eps'.format(k))
      fig.savefig('cctbx_sym_and_simple_vs_probe.png')
      plt.show()
      print 'Linear fit info: PROBE = {} * CCBTX + {}'.format(m,b)
      print ignore_str
      print '-'*50
    else:
      print 'Load data before attempting to plot'

def run():
  print 'Start'
  test_results = Collecting_clashscore_data()
  test_results.set_working_path()
  test_results.get_test_data()
  test_results.plot_reference(ignore_delta=10)
  test_results.plot_total(ignore_delta=10)
  print 'Done...'

if __name__=='__main__':
  run()
