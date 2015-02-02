from __future__ import division
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import cPickle as pickle
import numpy as np
import os, sys

try:
  from scipy.stats import linregress
except ImportError:
  print "Didn't load linregress from scipy.stats"


class Collecting_clashscore_data(object):

  def __init__(self):
    """
    Collect and analyse clashscore test data

    good test result should have only one line
      x1::x2::x3::x4::x5
      Where:
        x1 (0): pdb_file_name
        x2 (1): macro_molecule_cctbx_clashscore
        x3 (2): symmetry_cctbx_clashscore
        x4 (3): total_cctbx_clashscore
        x5 (4): probe_clashscore

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
    self.queue_data_path = "/net/cci/youval/work/work/Clashes/queue_clash_compare"
    # containers for results
    self.data_sym_clashes = {}
    self.data_outliers = {}     # when cctbx very different than probe
    #
    self.structure_where_pdb_not_found = []
    self.structure_with_error_when_processing_phenix_clashscore = []
    self.other_issues_with_results = []

  def set_working_path(self,all_data=True):
    """
    Set working and data folders
    When all_data=False use the data for the test on macro molecule (not
    the complete model)
    C:\Phenix\Dev\Work\work\Clashes\CCTBX_PROBE_compare
    """
    r1 = r'C:\Users\Youval\Google Drive\Documents\LBNL\phenix\publications'
    r1 += r'\news letter\clash score\related code'
    data_1 = r1 + r'\pdb_scan_result_macro_molecule_files'
    data_2 = r1 + r'\pdb_scan_result_complete_files'
    osType = sys.platform
    if osType.startswith('win'):
      assert os.path.isdir(data_1)
      assert os.path.isdir(data_2)
      if all_data:
        self.working_path = data_2
        print data_2
      else:
        self.working_path = data_1
        print data_1
    else:
      path = '/net/cci/youval/work/work/Clashes'
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

    "data_dict" is a dictionary for the "clean_data"
    "data" contains all data collected
    """
    have_data = os.path.isfile(self.data_file_name)
    have_data &= os.path.isfile(self.data_dict_file_name)
    if have_data:
      print 'using existing data files'
      self.data = pickle.load(open(self.data_file_name,'r'))
      self.clean_data = pickle.load(open(self.clean_data_file_name,'r'))
      self.clean_data = [x for x in self.clean_data if x[1] >= 0 ]
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
          if d[1] >= 0:
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
    for d in self.data:
      pdb_id = d[0]
      if d[1] == -1:
        self.structure_where_pdb_not_found.append(pdb_id)
      elif d[1] == -2:
        self.structure_with_error_when_processing_phenix_clashscore.append(pdb_id)
      elif d[1] == -3:
        self.other_issues_with_results.append(pdb_id)
    print "structure_where_pdb_not_found: ",len(self.structure_where_pdb_not_found)
    print "structure_with_error_when_processing_phenix_clashscore: ",\
      len(self.structure_with_error_when_processing_phenix_clashscore)
    n_to_print = min(10,len(self.structure_with_error_when_processing_phenix_clashscore))
    print self.structure_with_error_when_processing_phenix_clashscore[:n_to_print]
    print "other_issues_with_results: ",len(self.other_issues_with_results)
    print self.other_issues_with_results
    print '-'*50

  def plot_reference(self,ignore_delta=100):
    """
    Compare CCBTX macro molecule non-bonded clashscore to PROBE clashscore

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
      cctbx_prob = [(x[4],x[1]) for x in self.clean_data
                    if (abs(x[1]-x[4])<ignore_delta)]
      outliers = [x[0] for x in self.clean_data
                  if (abs(x[1]-x[4])>20) and (x[4] < 20)]

      n_ignored = len(self.clean_data)-len(cctbx_prob)
      n_all = len(self.clean_data)
      print '(macro. mol.) Number of data points ignored: ',n_ignored
      print outliers
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
      plt.xticks(fontsize=fontsize - 2)
      plt.yticks(fontsize=fontsize - 2)
      # plt.title(
      #   'CCBTX macro molecule vs PROBE non-bonded clashscore',
      #   fontsize=fontsize)
      plt.ylabel('PROBE clashscore',fontsize=fontsize)
      plt.xlabel('Non-bonded overlaps per 1000 atoms',fontsize=fontsize)
      ignore_str = 'Ignore {} of {} outlier points where\n'
      ignore_str +=  'abs(macro_mol cctbx_score - probe_score) >= {}'
      ignore_str = ignore_str.format(n_ignored,n_all,ignore_delta)
      # these are matplotlib.patch.Patch properties
      props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
      # plt.text(2,max_y - 2,ignore_str,
      #          fontsize=fontsize,verticalalignment='top',bbox=props)
      plt.xlim([0,max_x])
      plt.ylim([0,max_y])
      plt.show()
      fig.savefig('cctbx_macro_mol_vs_probe.png')
      print 'Linear fit info: PROBE = {} * CCBTX + {}'.format(m,b)
      print ignore_str
      print '-'*50
    else:
      print 'Load data before attempting to plot'

  def plot_total(self,ignore_delta=100):
    """
    Compare CCBTX total_cctbx_clashscore to PROBE clashscore

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
      assert len(self.clean_data[0]) == 5
      cctbx_prob = [(x[4],x[3]) for x in self.clean_data
                    if abs(x[4]-x[3])<ignore_delta]
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
      plt.xticks(fontsize=fontsize - 2)
      plt.yticks(fontsize=fontsize - 2)
      # plt.title(
      #   'CCBTX total vs PROBE non-bonded clashscore',fontsize=fontsize)
      plt.ylabel('PROBE clashscore',fontsize=fontsize)
      plt.xlabel( 'Non-bonded overlaps per 1000 atoms', fontsize=fontsize)
      ignore_str = 'Ignore {} of {} outlier points where\n'
      ignore_str +=  'abs(clashscore difference) >= {}'
      ignore_str = ignore_str.format(n_ignored,n_all,ignore_delta)
      # these are matplotlib.patch.Patch properties
      props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
      # plt.text(2,max_y - 2,ignore_str,
      #          fontsize=fontsize,verticalalignment='top',bbox=props)
      plt.xlim([0,max_x])
      plt.ylim([0,max_y])
      plt.show()
      fig.savefig('cctbx_total_vs_probe.png')
      print 'Linear fit all data, info: PROBE = {} * CCBTX + {}'.format(m,b)
      print ignore_str
      print '-'*50
    else:
      print 'Load data before attempting to plot'

  def plot_sym(self,ignore_delta=100):
    """
    Plot plot of CCBTX all clashscore and symmetry clashscore vs. PROBE
    clashscore

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
      assert len(self.clean_data[0]) == 5
      # x[4] : probe
      # x[3] : cctbx_total_vs_probe
      # x[2] : symmetry_cctbx_clashscore
      cctbx_prob = [(x[4],x[3],x[2]) for x in self.clean_data
                    if abs(x[4]-x[3])<ignore_delta]
      n_ignored = len(self.clean_data)-len(cctbx_prob)
      n_all = len(self.clean_data)
      print 'Number of data points ignored: ',n_ignored
      cctbx_prob.sort()
      # cctbx
      cctbx_score = [x[1] for x in cctbx_prob]
      cctbx_sym = [x[2] for x in cctbx_prob]
      probe_score = [x[0] for x in cctbx_prob]
      # Get linear fitting parameters
      x_fit = [0,max_x]
      # m,b = plb.polyfit(cctbx_score, probe_score, 1)
      m,b ,r_value,_,_ = linregress(cctbx_score, probe_score)
      print 'r-value: ',r_value
      y_fit = [b,m * max_x + b]
      print 'Linear fit info: PROBE = {} * CCBTX + {}'.format(m,b)
      print '-'*50
      # create plot
      # gr = 1.61803398875 	# Golden ratio
      plt.close('all')
      gr = 1
      h = 4				# figure height
      w = gr*2*h			# figure width
      # setup subplots
      fig = plt.figure(figsize=(8.3,8.2))
      # gs = gridspec.GridSpec(2,1,width_ratios=[1,1],height_ratios=[2,1])
      gs = gridspec.GridSpec(2,1,height_ratios=[2,1])
      # gs.update(left=0.05, right=0.48, wspace=0.05)
      ax1 = plt.subplot(gs[0,0])
      ax2 = plt.subplot(gs[1,0])
      ax1.plot(cctbx_score,probe_score,'.b',x_fit,y_fit,'y',linewidth=2)
      ax1.tick_params(axis='both',labelsize=fontsize)
      ax1.ticklabel_format(axis='both',labelsize=fontsize)
      # ax1.set_title(
      #   'Clashscores and symmetry related clashes',
      #   fontsize=fontsize)
      ax1.set_ylabel('PROBE clashscore',fontsize=fontsize)
      ax1.set_xlabel('Non-bonded overlaps per 1000 atoms',fontsize=fontsize)
      # these are matplotlib.patch.Patch properties
      props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
      text_params = {'fontsize':fontsize,'verticalalignment':'top','bbox':props}
      ax1.text(2,max_y - 2,'a',**text_params)
      ax1.set_xlim([0,max_x])
      ax1.set_ylim([0,max_y])

      # second plot
      n_files = np.arange(0,len(cctbx_sym))
      ax2.plot(cctbx_sym,n_files,'.g')
      ax2.text(2,max_y - 2,'b',**text_params)
      ax2.set_ylabel('Files',fontsize=fontsize)
      ax2.set_xlabel('CCTBX symmetry clashscore',fontsize=fontsize)
      ax2.tick_params(axis='both',labelsize=fontsize)
      ax2.set_ylim([0,len(cctbx_sym)])
      ax2.set_xlim([0,max_y])
      ax2.set_yticks((5000,20000))
      ax2.set_yticklabels(('5k','20k'))
      plt.show()
      # fig.savefig('cctbx_simple_vs_probe.png')
    else:
      print 'Load data before attempting to plot'

  def hist_sym(self,prob_limit=1000):
    """
    CCTBX symmetry clashscore histogram

    - Considering clashes where PROBE clashscore is <= prob_limit
    - All symmetry clashscores larger than 9 will be added to the last bin

    Args:
      ignore_delta (float): ignore outlier points where
        abs(cctbx_score - probe_score) > ignore_delta
    """
    prefix = 'probe_{}_'.format(prob_limit)
    # n_bins = 50
    n_bins = range(0,11,1)
    # n_bins = range(0,ignore_delta,1)
    fontsize = 20
    fig = plt.figure(figsize=(12,8.3))
    print dir(fig)

    # get data
    assert len(self.clean_data[0]) == 5
    # x[2] : symmetry_cctbx_clashscore
    # x[4] : MolProbity (PROBE) clashscore
    cctbx_sym = [x[2] for x in self.clean_data if x[4]<=prob_limit]
    t = lambda x: x*(x<=9) + 9.5*(x>9)
    cctbx_sym = [t(x) for x in cctbx_sym]
    # t = lambda x,d,l: (x[2]<d) and (x[4]<l)
    # cctbx_sym = [x[2] for x in self.clean_data if t(x,ignore_delta,prob_limit)]
    # cctbx_sym = [x[2] for x in self.clean_data if (x[2] <= ignore_delta)]
    n_ignored = len(self.clean_data)-len(cctbx_sym)
    print 'Histogram - Number of data points ignored: ',n_ignored
    print 'Histogram - total number of points: ',len(cctbx_sym)
    # s = 'Histogram - Number sym clashscore larger then {}: {}'
    # high_sym_clash = [x for x in self.clean_data if x[2] > ignore_delta]
    # print s.format(ignore_delta,len(high_sym_clash))


    # histogram our data with numpy
    n, bins, patches = plt.hist(cctbx_sym, bins=n_bins)
    # we need to normalize the data to 0..1 for the full
    # range of the colormap
    import matplotlib.colors as colors
    import matplotlib.cm as cm
    fracs = np.log(n.astype(float))/n.max()
    norm = colors.Normalize(fracs.min(), fracs.max())

    for thisfrac, thispatch in zip(fracs, patches):
        color = cm.jet(norm(thisfrac))
        thispatch.set_facecolor(color)

    # plt.title('CCTBX symmetry related clashscore',fontsize=fontsize)
    plt.xlabel('Non-bonded symmetry overlaps per 1000 atoms',fontsize=fontsize)
    plt.ylabel('Number of PDB structures',fontsize=fontsize)
    plt.xlim((0,max(0.1,bins.max())))
    plt.xticks([3,6,9,10],['3','6','9','$\infty$'],fontsize=fontsize - 2)
    plt.yticks(fontsize=fontsize - 2)
    plt.yscale('log')
    plt.show()
    fig.savefig(prefix + 'cctbx_sym_clashscore_hist.png',fontsize=fontsize)

  def get_sym_clashscore(self,larger_then = 10):
    """ Collect all structures with clashscore due to symmetry is
    larger than 'larger_then' """
    self.data_sym_clashes = {x[0]:x for x in self.clean_data
                             if x[2] > larger_then}
    msg = 'Structures where the clashscore due to symmetry is larger than {}: {}'
    print msg.format(larger_then,len(self.data_sym_clashes))
    show_table(n=10,dict_to_print=self.data_sym_clashes)


  def get_outliers(self,delta = 10,min_clashscore=20):
    """ Collect all structures with clashscore where the difference between
    the CCTBX and the PROBE clashscores is 'larger then delta' """
    mc = min_clashscore
    self.data_outliers = {x[0]:x for x in self.clean_data
                          if (abs(x[1]-x[4]) > delta)
                          and
                          ((x[1] <= mc) or (x[4] <= mc))}

    msg = 'Structures where CCTBX and PROBE clashscore difference'
    msg += ' is larger than {} and one of them is smaller than {}: {}'
    print msg.format(delta,min_clashscore,len(self.data_outliers))
    show_table(n=10,dict_to_print=self.data_outliers,rows_to_print=10)

def show_table(n,dict_to_print=None,list_to_print=None,rows_to_print=0):
  """ Print as table with 'n' items in each raw """
  l = len(dict_to_print)
  m = l // n
  dn = n - (l - n*m)
  if dict_to_print:
    data = dict_to_print.keys()
  else:
    assert list_to_print
    data = list_to_print
  data.extend(['',] *dn)
  if rows_to_print == 0: rows_to_print = 100
  rows_to_print = min(rows_to_print,len(data)//n)
  for i in range(rows_to_print):
    print data[i*n:(i+1)*n]
  print '.'*100

def run():
  print 'Start'
  test_results = Collecting_clashscore_data()
  test_results.set_working_path(all_data=False)
  # test_results.set_working_path(all_data=True)
  test_results.get_test_data()
  test_results.plot_reference()
  # test_results.plot_sym(ignore_delta=100)
  # for fig. 4
  # test_results.hist_sym(prob_limit=10000)
  # for fig. 5
  test_results.hist_sym(prob_limit=0)
  # test_results.plot_total(ignore_delta=100)
  test_results.get_sym_clashscore(larger_then = 10)
  test_results.get_outliers(delta=10,min_clashscore=0)
  print 'Done...'

if __name__=='__main__':
  run()
