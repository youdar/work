from __future__ import division
import matplotlib.pyplot as plt
from pprint import pprint
import cPickle as pickle
import numpy as np
import os
import sys
import re

"""
File data input example:
#  PDB code |Reported in PDB  | Calc from NCS   |  ASU initial    |   ASU final     | Res.   |   NCS   |  Solvent |    Data      | Year   | Use  | Geo. | Trans. | time (sec)
#           | r-work | r-free | r-work | r-free | r-work | r-free | r-work | r-free |        | copies  | fraction | completeness |        | NCS  | Rest | refine |
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
     tttt   | -1.00  | -1.00  |  0.57  |  0.58  |  0.33  |  0.34  |  0.00  |  0.00  | -1.00  |    3    |   0.95   |     1.00     |  None  | True |False |  True  |   14
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Rotation,Translation difference      |  Rxx   |  Rxy   |  Rxz   |  Ryx   |  Ryy   |  Ryz   |  Rzx   |  Rzy   |  Rzz   |   Tx   |   Ty   |   Tz
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  Transform number: 1                  | 0.0454 | 0.0219 |-0.0169 |-0.0278 | 0.0410 |-0.0222 |-0.0101 |-0.0036 | 0.0030 | 0.0006 | 0.0005 |-0.0004
  Transform number: 2                  |-0.0439 |-0.0115 |-0.0184 |-0.0579 |-0.0646 |-0.0046 | 0.0465 |-0.0623 | 0.0493 | 0.0003 |-0.0010 | 0.0001

=============================================================================================================================================================================

Posiotion in the data list:
      0         1	     2	      3	        4	     5	      6        7       8         9       10        11            12         13       14     15     16         17


Data to excel output
--------------------
   0 (0)|  1 (10)        | 2(13)|   3(9)     |      4(12)       |   5(11)          |
pdb code| num ncs copies | year | resolution | data completeness| solvent fraction |


 6(1)             |     7(2)          |     8(7)(14)            |    9(8)(14)
R-work pdb header | R-free pdb header | R-work using strict ncs | R-free using strict ncs |

  10(7)(14)               |        11(8)(14)          |       12(7)(16)       |
R-work without strict ncs | R-free without strict ncs | R-work with transform |

  13(8)(16)           |    14(5)        |    15(6)        |   16(15)
R-free with transform | r work asu init | r free asu init |	use geometry restraints |

  17(16)      | 18(14)(14)| 19(14)(14)  |    20(14)(16)  |
use transform | time ncs  | time no ncs | time_transform |

  21
Compare to reported

"""

class results_collection(object):

  def __init__(self):
    self.data_records = []
    self.data_records_long = []
    self.data_records_strict_ncs_dict = {}
    self.data_records_without_strict_ncs_dict = {}
    self.data_records_dict = {}
    self.cols_inp = [
      'pdb_code',
      'r_work_pdb_reported',
      'r_free_pdb_reported',
      'r_work_single_ncs',
      'r_free_single_ncs',
      'r_work_asu_init',
      'r_free_asu_init',
      'r_work_asu_final',
      'r_free_asu_final',
      'resolution',
      'num_ncs_copies',
      'solvent_fraction',
      'data_completeness',
      'year',
      'use_strict_ncs',
      'use_geometry_restraints',
      'use_transforms',
      'time']

    self.cols_names = [
      'pdb_code',
      'num_ncs_copies',
      'year',
      'resolution',
      'data_completeness',
      'solvent_fraction',
      'r_work_pdb_header',
      'r_free_pdb_header',
      'r_work_asu_ncs',
      'r_free_asu_ncs',
      'r_work_asu_no_ncs',
      'r_free_asu_no_ncs',
      'r_work_asu_transform',
      'r_free_asu_transform',
      'r_work_asu_init',
      'r_free_asu_init',
      'use_geometry_restraints',
      'use_transforms',
      'time_ncs',
      'time_no_ncs',
      'time_transform']

    self.cols_names_for_table1 = [x.replace('_',' ') for x in self.cols_inp]
    self.cols_names_for_table2 = [x.replace('_',' ') for x in self.cols_names]
    self.map_to_ncs = {0:0,1:6,2:7,5:14,6:15,13:2,9:3,10:1,11:5,12:4,15:16,
                       16:17,7:8,8:9,14:18}
    self.map_to_no_ncs = {0:0,1:6,2:7,5:14,6:15,13:2,9:3,10:1,11:5,12:4,15:16,
                          16:17,7:10,8:11,14:19}
    self.map_to_transform = {0:0,1:6,2:7,5:14,6:15,13:2,9:3,10:1,11:5,12:4,
                             15:16,16:17,7:12,8:13,14:20}
    self.float_type_records = [1,2,3,4,5,6,7,8,9,11,12,17]
    self.int_type_records = [10,13]
    self.bool_type_records = [14,15,16]
    self.str_type_records = [0]

    class sort_inp():
      def __init__(self):
        """ Parameter number to sort by """
    sort_inp.__dict__.update(zip(self.cols_inp, range(len(self.cols_inp))))
    self.sort_inp = sort_inp

  def read_filenames(self):
    use_rec = self.sort_inp
    # Read files in current directory
    files = os.listdir(os.getcwd())
    # collect only non-empty files that starts with log_
    files = [x for x in files if x.startswith('log_') and os.stat(x).st_size>0]
    print 'Number of log files with data: {}'.format(len(files))
    # prepare regular expression search of all lines that start with pdb code
    regex = re.compile('     ....   \|')
    for fn in files:
      pdb_code = fn[4:]
      d = open(fn, "r").readlines()
      for ln in d:
        # print info and warnings
        if ln.startswith('Warning'):
          print '{0}: {1}'.format(pdb_code,ln.strip())
        if ln.startswith('Using pdb file from local machine'):
          print 'Used {} from local machine'.format(pdb_code)
        if ln.startswith('Traceback'):
          print '{0}: {1}'.format(pdb_code,'Problem processing')
        # process files with good data
        if re.match(regex,ln[:13]):
          data = [x.strip() for x in ln.split('|')]
          # convert data values to proper type
          for indx in self.float_type_records:
            data[indx] = float(data[indx])
          for indx in self.int_type_records:
            if not data[indx]: data[indx] = 0
            else: data[indx] = int(data[indx])
          for indx in self.bool_type_records:
            data[indx] = (data[indx] == 'True')
          # collect data as is
          self.data_records.append(data)
          # combine all data for each pdb file
          if self.data_records_dict.has_key(data[0]):
            existing_record = self.data_records_dict[data[0]]
          else:
            existing_record = len(self.cols_names)*[None,]
          if data[use_rec.use_strict_ncs] and not data[use_rec.use_transforms]:
            for map_from,map_to in self.map_to_ncs.iteritems():
              if not existing_record[map_to]:
                existing_record[map_to] = data[map_from]
          if not data[use_rec.use_strict_ncs] and not data[use_rec.use_transforms]:
            for map_from,map_to in self.map_to_no_ncs.iteritems():
              if not existing_record[map_to]:
                existing_record[map_to] = data[map_from]
          if data[use_rec.use_transforms]:
            for map_from,map_to in self.map_to_transform.iteritems():
              if not existing_record[map_to]:
                existing_record[map_to] = data[map_from]
          self.data_records_dict[data[0]] = existing_record
          # build other dictionaries
          if data[use_rec.use_strict_ncs]:
            self.data_records_strict_ncs_dict[data[0]] = data
          else:
            self.data_records_without_strict_ncs_dict[data[0]] = data
    # create a list version of the dictionary
    for key,val in self.data_records_dict.iteritems():
      self.data_records_long.append(val)
    print 'There are {} good data records'.format(len(self.data_records))

    # organize data
    # Sort records
    # self.data_records.sort(key=lambda x:x[use_rec.time])
    print 'done with collection'
    pickle.dump(self.data_records,open('data_records_v2','w'))
    pickle.dump(self.data_records_long,open('data_records_long_v2','w'))
    pickle.dump(self.data_records_dict,open('data_records_dict_v2','w'))
    pickle.dump(
      self.data_records_strict_ncs_dict,open('data_records_strict_ncs_dict_v2','w'))
    pickle.dump(
      self.data_records_without_strict_ncs_dict,
      open('data_records_without_strict_ncs_dict_v2','w'))

  def get_data_from_files(self):
    """  Read data from files    """
    self.data_records = pickle.load(open('data_records_v2','r'))
    self.data_records_long = pickle.load(open('data_records_long_v2','r'))
    self.data_records_dict = pickle.load(open('data_records_dict_v2','r'))
    self.data_records_strict_ncs_dict = pickle.load(
      open('data_records_strict_ncs_dict_v2','r'))
    self.data_records_without_strict_ncs_dict = pickle.load(
      open('data_records_without_strict_ncs_dict_v2','r'))
    print 'Got the data from files'

  def plot_results_1(self):
    """
    Plot a result summery:
    r-work (refined using strict-ncs) vs. r-work (refined without strict-ncs)

    xerr , yerr are the difference between r-work and r-free for the two
    refinement methods

    The size of the circles indicate the size of (xerr - yerr). It is blue
    when the strict_ncs refinement make the (r_work - r_free) smaller.
    """
    sort_by = self.sort_inp
    # Collect data points that have results with and without strict_ncs
    x_ncs = []; y_no_ncs =[]
    x2 = []; y2 = []
    time_ncs = []; time_no_ncs = []
    pdb_code = []
    for x in self.data_records_strict_ncs_dict:
      if self.data_records_without_strict_ncs_dict.has_key(x):
        x_ncs.append(self.data_records_strict_ncs_dict[x]
                           [sort_by.r_work_asu_final])
        y_no_ncs.append(self.data_records_without_strict_ncs_dict[x]
                                  [sort_by.r_work_asu_final])
        x2.append(self.data_records_strict_ncs_dict[x]
                      [sort_by.r_free_asu_final])
        y2.append(self.data_records_without_strict_ncs_dict[x]
                      [sort_by.r_free_asu_final])
        time_ncs.append(self.data_records_strict_ncs_dict[x]
                      [sort_by.time])
        time_no_ncs.append(self.data_records_without_strict_ncs_dict[x]
                      [sort_by.time])
        pdb_code.append(x)

    # variable error bar values
    xerr = [x-y for (x,y) in zip(x_ncs,x2)]
    yerr = [x-y for (x,y) in zip(y_no_ncs,y2)]

    # get larges r-values for the 45 degrees line
    maxval = max(x_ncs + y_no_ncs) *1.05

    # Convert lists to numpy arrays
    x_ncs = np.array(x_ncs)
    y_no_ncs = np.array(y_no_ncs)
    xerr = np.array(xerr)
    yerr = np.array(yerr)
    # First illustrate basic pyplot interface, using defaults where possible.
    plt.figure()
    colors = []
    delta_err = [(x-y) for (x,y) in zip (x_ncs,y_no_ncs)]
    s = [5+1000*abs(x) for x in delta_err]
    for i,x in enumerate(delta_err):
      if x > 0:
        colors.append('b')
      else:
        colors.append('y')
        print 'Better r-free without NCS: ',pdb_code[i]
    # add points for size reference : Difference of 0, 0.05, 0.1
    s.extend([105,55,5])
    colors.extend(['g','g','g'])
    x_ncs = np.append(x_ncs,[.02,.02,.02])
    d = maxval - 0.55
    y_no_ncs = np.append(y_no_ncs,[.4 + d,.45 + d,.5 + d])
    plt.scatter(x_ncs,y_no_ncs,s,c=colors)
    plt.plot([0,maxval],[0,maxval])
    plt.xlabel('r-work Refinement using strict-ncs')
    plt.ylabel('r-work Refinement without strict-ncs')
    plt.title('Compare r-work values and stric-ncs influence on r-free')
    # plot reference point, to indicate the meaning of size
    plt.text(0.035,0.49 + d, 'The same r_work - r_free value',fontsize=14)
    plt.text(0.035,0.44 + d, '0.05 difference',fontsize=14)
    plt.text(0.035,0.39 + d, '0.1 difference',fontsize=14)

    plt.xlim(0,maxval)
    plt.ylim(0,maxval)
    plt.show()

  def plot_results_2(self):
    """
    Plot a result summery:
    time for refinement with strict ncs vs. without
    """
     # Collect data points that have results with and without strict_ncs
    time_ncs = []; time_no_ncs = []
    pdb_code = []
    for x in self.data_records_strict_ncs_dict:
      if self.data_records_without_strict_ncs_dict.has_key(x):
        time_ncs.append(self.data_records_strict_ncs_dict[x]
                      [self.sort_inp.time])
        time_no_ncs.append(self.data_records_without_strict_ncs_dict[x]
                      [self.sort_inp.time])
        pdb_code.append(x)

    # print outliers
    d = 0.1
    for (x,y,code) in zip(time_ncs,time_no_ncs,pdb_code):
      if abs(x-y)/x > d:
        print 'Time difference is {0:.0f}% for {1}'.format(100*(x-y)/x,code)
    # get larges r-values for the 45 degrees line
    maxval = max(time_ncs + time_no_ncs) *1.05

    # Convert lists to numpy arrays
    x_ncs = np.array(time_ncs)
    y_no_ncs = np.array(time_no_ncs)
    # plot
    plt.figure()
    plt.plot(x_ncs,y_no_ncs,'o',[0,maxval],[0,maxval])
    plt.xlabel('Time[sec] Refinement using strict-ncs')
    plt.ylabel('Time[sec] Refinement without strict-ncs')
    plt.title('Looking at refinement time when using stric_ncs')
    plt.xlim(0,maxval)
    plt.ylim(0,maxval)
    plt.show()

  def save_csv_table_to_file(self):
    """
    Save a the collected data in a csv file:

    The data in the file is one line per test with the following records order:
    'pdb_code','r_work_pdb_reported','r_free_pdb_reported',
    'r_work_single_ncs','r_free_single_ncs','r_work_asu_init',
    'r_free_asu_init','r_work_asu_final','r_free_asu_final','resolution',
    'num_ncs_copies','solvent_fraction','data_completeness','year',
    'use_strict_ncs','use_geometry_restraints','use_transforms','time'
    """
    file_name = 'ncs_refinement_results_v2.csv'
    f = open(file_name,'w')
    f.write(','.join(self.cols_names_for_table1))
    f.write('\n')
    for rec in self.data_records:
      rec = [str(x) for x in rec]
      rec += '\n'
      f.write(','.join(rec))
    f.close()
    print 'data was saved to: ',file_name

  def save_csv_table2_to_file(self):
    """
    Save a the collected data in a csv file:

    The data in the file is one line per test with the following records order:

    'pdb_code','num_ncs_copies','year','resolution','data_completeness',
    'solvent_fraction','r_work_pdb_header','r_free_pdb_header',
    'r_work_asu_ncs','r_free_asu_ncs','r_work_asu_no_ncs','r_free_asu_no_ncs',
    'r_work_asu_transform','r_free_asu_transform','r_work_asu_init',
    'r_free_asu_init','use_geometry_restraints','use_transforms',
    'time_ncs','time_no_ncs','time_transform'
    """
    file_name = 'ncs_refinement_results_table2_v2.csv'
    table_title = self.cols_names_for_table2
    f = open(file_name,'w')
    f.write(','.join(self.cols_names_for_table2))
    f.write('\n')
    for rec in self.data_records_long:
      rec = [str(x) for x in rec]
      rec += '\n'
      f.write(','.join(rec))
    f.close()
    print 'data was saved to: ',file_name


  def get_list_of_unprocessed_files(self):
    """
    Check which of the PDB files, from our initial list, are not included in
    the results
    """
    pdb_code_set = {
      '3dar', '1vcr', '1r2j', '1a37', '1llc', '1tnv', '1tdi', '1w39', '1ny7',
      '1ddl', '1c8n', '2bfu', '4gmp', '3vbr', '3vbu', '3vbo', '4jgy', '3es5',
      '3nop', '3not', '3nou', '3bcc', '1bcc', '1z7s', '6msf', '2iz8', '7msf',
      '2izn', '2c50', '2c51', '2iz9', '2c4y', '2c4z', '5msf', '2c4q', '2bu1',
      '3raa', '3oah', '3ra2', '3ra9', '3ra8', '3ra4', '3qpr', '1ei7', '1a34',
      '3chx', '2wbh', '2fz1', '2fz2', '2gh8', '1wcd', '3fbm', '4gb3', '1laj',
      '3vbh', '1dzl', '3hag', '4iv3', '1js9', '3n7x', '4gh4', '4jgz', '3tn9',
      '4iv1', '1vb2', '1vb4', '1vak', '3s4g', '2buk', '1x36', '4bcu', '1b35',
      '2wzr', '1k5m', '2bq5', '1zba', '1pgw', '3vbs', '1x35', '3vbf', '1pgl',
      '4fsj', '4fte', '4fts', '2e0z', '4ftb', '2w4y', '2w4z', '2qzv', '3vdd',
      '3p0s', '1qjx', '1qjy', '1qju', '3r0r', '2bs1', '2ztn', '1x9t', '2zzq',
      '1x9p', '4aqq', '1za7', '4ar2', '2wws', '2xpj', '4hl8', '3ntt', '2vf1',
      '3ux1', '2xgk', '2izw', '3cji', '4gbt', '2vq0', '4g93', '2g34', '2qij',
      '2g33', '1f2n', '4g0r', '1ng0', '2ws9', '2xbo', '2wff', '1wce', '1dwn',
      '2vf9', '3zfe', '3zff', '3zfg', '2x5i', '1h8t', '3lob', '4ang', '2gtl',
      '2qqp', '1f8v', '1m1c', '1lp3', '4aed', '3e8k', '1uf2', '1ohg', '1ohf',
      '3s6p', '3kz4', '4f5x', '1vsz'}

    without_strict_ncs_set=set(self.data_records_without_strict_ncs_dict.keys())
    with_strict_ncs_set=set(self.data_records_strict_ncs_dict.keys())
    print ''
    print '  files with some records missing'
    print '---------------------------------'
    outdata1 = pdb_code_set - without_strict_ncs_set
    outdata2 = pdb_code_set - with_strict_ncs_set
    outdata = list(outdata1.union(outdata2))
    n = 5
    l = len(outdata)
    for i in range(0,l,n):
      e = min(i+n,l)
      s = ['{}']*(e-i)
      s = ', '.join(s)
      print s.format(*outdata[i:e])

if __name__=='__main__':
  # path_to_log_files = "/net/cci/youval/Work/work/NCS/junk/pdb_test/queue_job"
  path_to_log_files = r"C:\Phenix\Dev\Work\work\NCS\junk\pdb_test\queue_job"
  current_path = os.getcwd()
  os.chdir(path_to_log_files)
  process_results = results_collection()
  # need to read_filenames only for new runs, it overwrites existing data files
  process_results.read_filenames()
  # process existing files
  process_results.get_data_from_files()
  process_results.save_csv_table_to_file()
  process_results.save_csv_table2_to_file()
  # plot
  process_results.plot_results_1()
  process_results.plot_results_2()
  # check which files are not processed yet
  process_results.get_list_of_unprocessed_files()
  os.chdir(current_path)
  print 'Done...'
