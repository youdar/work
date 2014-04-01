from __future__ import division
import matplotlib.pyplot as plt
import cPickle as pickle
import numpy as np
import os
import sys
import re

"""
Output example:
#  PDB code |Reported in PDB  | Calc from NCS   |  ASU initial    |   ASU final     | Res.   | Year   | Use NCS  |Geo. Rest. | time (sec)
#           | r-work | r-free | r-work | r-free | r-work | r-free | r-work | r-free |        |        |          |           |
#----------------------------------------------------------------------------------------------------------------------------------------
     1vcr   |  0.38  |  0.35  |  0.44  |  0.45  |  0.43  |  0.41  |  0.82  |  0.87  |  9.50  |  2004  |  False   |   True    |   136

Posiotion in the data list:
      0         1	  2	   3	   4	     5	      6        7       8         9       10        11         12          13
"""

class results_collection(object):

  def __init__(self):
    self.data_records = []
    self.data_records_strict_ncs_dict = {}
    self.data_records_dict = {}
    self.cols_names = [
      'pdb_code','r_work_pdb_reported','r_free_pdb_reported',
      'r_work_pdb_ncs','r_free_pdb_ncs',
      'r_work_asu_init','r_free_asu_init',
      'r_work_asu_final','r_free_asu_final',
      'resolution','year','use_strict_ncs',
      'use_geometry_restraints','time']
    self.cols_names_for_table = [
      'pdb code','r-work pdb reported',
      'r-free pdb reported',
      'r-work pdb ncs','r-free pdb ncs',
      'r-work asu init','r-free asu init',
      'r-work asu final','r-free asu final',
      'resolution','year ','use strict ncs',
      'use geometry restraints','time']


  def read_filenames(self):
    class sort_by():
      # def __init__(self):
        """ Parameter number to sort by"""
    sort_by.__dict__.update(zip(self.cols_names, range(len(self.cols_names))))
    assert sort_by.time == 13
    # Read files in current directory
    files = os.listdir(os.getcwd())
    # collect only non-empty files that starts with log_
    files = [x for x in files if x.startswith('log_') and os.stat(x).st_size>0]
    print 'Number of log files with data: {}'.format(len(files))
    # prepare regular expression serch
    regex = re.compile('     ....   \|')
    for fn in files:
      pdb_code = fn[4:]
      d = open(fn, "r").readlines()
      for ln in d:
        if ln.startswith('Warning'):
          print '{0}: {1}'.format(pdb_code,ln.strip())
        if ln.startswith('Traceback'):
          print '{0}: {1}'.format(pdb_code,'Problem processing')
        if re.match(regex,ln[:13]):
          data = [x.strip() for x in ln.split('|')]
          for indx in [1,2,3,4,5,6,7,8,9]:
            # convert r-values to float
            data[indx] = float(data[indx])
          # convert year and time to integers
          data[10] = round(float(data[10]))
          data[13] = round(float(data[13]))
          # print data
          self.data_records.append(data)
          if data[sort_by.use_strict_ncs] == 'True':
            self.data_records_strict_ncs_dict[data[0]] = data
          else:
            self.data_records_dict[data[0]] = data
    print 'There are {} good data records'.format(len(self.data_records))
    # Sort records
    # self.data_records.sort(key=lambda x:x[sort_by.time])
    print 'done with collection'


  def plot_results_1(self):
    """
    Plot a result summery:
    (error bars reference:
    http://matplotlib.org/1.2.1/examples/pylab_examples/errorbar_demo.html)

    Reported r-work vs. (r-work of ASU,  with (blue) and without (red)
    strict_ncs)

    plot error bars that represent (r_work - r_free). Horizontal error bar for
    published values, vertical for values after refinement.
    """
    # example data
    x = np.arange(0.1, 4, 0.5)
    y = np.exp(-x)

    # example variable error bar values
    yerr = 0.1 + 0.2*np.sqrt(x)
    xerr = 0.1 + yerr

    # First illustrate basic pyplot interface, using defaults where possible.
    plt.figure()
    plt.errorbar(x, y, xerr=0.2, yerr=0.4)
    plt.title("Simplest errorbars, 0.2 in x, 0.4 in y")

    # Now switch to a more OO interface to exercise more features.
    fig, axs = plt.subplots(nrows=2, ncols=2, sharex=True)
    ax = axs[0,0]
    ax.errorbar(x, y, yerr=yerr, fmt='o')
    ax.set_title('Vert. symmetric')

    # With 4 subplots, reduce the number of axis ticks to avoid crowding.
    ax.locator_params(nbins=4)

    ax = axs[0,1]
    ax.errorbar(x, y, xerr=xerr, fmt='o')
    ax.set_title('Hor. symmetric')

    ax = axs[1,0]
    ax.errorbar(x, y, yerr=[yerr, 2*yerr], xerr=[xerr, 2*xerr], fmt='--o')
    ax.set_title('H, V asymmetric')

    ax = axs[1,1]
    ax.set_yscale('log')
    # Here we have to be careful to keep all y values positive:
    ylower = np.maximum(1e-2, y - yerr)
    yerr_lower = y - ylower

    ax.errorbar(x, y, yerr=[yerr_lower, 2*yerr], xerr=xerr,
                fmt='o', ecolor='g', capthick=2)
    ax.set_title('Mixed sym., log y')

    fig.suptitle('Variable errorbars')

    plt.show()

  def plot_results_2(self):
    """
    Plot a result summery:

    r-work of ASU,  with  vs. without strict_ncs

    plot error bars that represent (r_work - r_free). Horizontal/Vertical
    error bars correspond to refinement with and without strict_ncs
    """

  def save_csv_table_to_file(self):
    """
    Save a the collected data in a csv file:

    The data in the file is one line per test with the following records order:
    'pdb_code','r_work_pdb_reported','r_free_pdb_reported',
    'r_work_pdb_ncs','r_free_pdb_ncs','r_work_asu_init','r_free_asu_init',
    'r_work_asu_final','r_free_asu_final','resolution','year','use_strict_ncs',
    'usd_geometry_restaints','time'
    """
    file_name = 'ncs_refinement_results.csv'
    f = open(file_name,'w')
    f.write(','.join(self.cols_names_for_table))
    f.write('\n')
    for rec in self.data_records:
      rec = [str(x) for x in rec]
      rec += '\n'
      f.write(','.join(rec))
    f.close()
    print 'data was saved to: ',file_name


if __name__=='__main__':
  # path_to_log_files = "/net/cci/youval/Work/work/NCS/junk/pdb_test/queue_job"
  path_to_log_files = r"C:\Phenix\Dev\Work\work\NCS\junk\pdb_test\queue_job"
  current_path = os.getcwd()
  os.chdir(path_to_log_files)
  process_results = results_collection()
  process_results.read_filenames()
  # process_results.save_csv_table_to_file()
  process_results.plot_results_1()
  os.chdir(current_path)
  print 'Done...'
