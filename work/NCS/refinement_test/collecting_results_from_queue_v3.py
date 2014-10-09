from __future__ import division
import matplotlib.pyplot as plt
import cPickle as pickle
import numpy as np
import os
import sys
import re

"""
File data input example:
#  PDB code |Reported in PDB  | Calc from NCS   |  ASU initial    |   ASU final     | Res.   |   NCS   |  Solvent |     Data     |  D/A ratio   |  D/A ratio   | Year   | Use  | Geo. | Trans. | time (sec)
#           | r-work | r-free | r-work | r-free | r-work | r-free | r-work | r-free |        | copies  | fraction | completeness |     ASU      |     NCS      |        | NCS  | Rest | refine |
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
     1vcr   | 0.3790 | 0.3530 | 0.4292 | 0.4486 | 0.4242 | 0.4382 | 0.3280 | 0.4407 |  9.50  |    5    |   0.93   |     1.00     |     0.21     |     1.06     |  2004  |False | True | False  |  1131

Using pdb file from local machine (not from MIRROR)
Using cif file from local machine (not from MIRROR)
#  PDB code |Reported in PDB  | Calc from NCS   |  ASU initial    |   ASU final     | Res.   |   NCS   |  Solvent |     Data     |  D/A ratio   |  D/A ratio   | Year   | Use  | Geo. | Trans. | time (sec)
#           | r-work | r-free | r-work | r-free | r-work | r-free | r-work | r-free |        | copies  | fraction | completeness |     ASU      |     NCS      |        | NCS  | Rest | refine |
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
     1vcr   | 0.3790 | 0.3530 | 0.4292 | 0.4486 | 0.4242 | 0.4382 | 0.3662 | 0.4119 |  9.50  |    5    |   0.93   |     1.00     |     0.21     |     1.06     |  2004  | True | True | False  |  928

Using pdb file from local machine (not from MIRROR)
Using cif file from local machine (not from MIRROR)
#  PDB code |Reported in PDB  | Calc from NCS   |  ASU initial    |   ASU final     | Res.   |   NCS   |  Solvent |     Data     |  D/A ratio   |  D/A ratio   | Year   | Use  | Geo. | Trans. | time (sec)
#           | r-work | r-free | r-work | r-free | r-work | r-free | r-work | r-free |        | copies  | fraction | completeness |     ASU      |     NCS      |        | NCS  | Rest | refine |
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
     1vcr   | 0.3790 | 0.3530 | 0.4292 | 0.4486 | 0.4242 | 0.4382 | 0.3558 | 0.4175 |  9.50  |    5    |   0.93   |     1.00     |     0.21     |     1.06     |  2004  | True | True |  True  |  1554
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Rotation,Translation difference      |  Rxx   |  Rxy   |  Rxz   |  Ryx   |  Ryy   |  Ryz   |  Rzx   |  Rzy   |  Rzz   |   Tx   |   Ty   |   Tz
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  Transform number: 1                  |-0.0132 | 0.0194 | 0.0100 |-0.0370 | 0.0708 | 0.0253 |-0.0481 |-0.0465 | 0.0563 | 0.8755 | 1.6418 | 3.8338
  Transform number: 2                  | 0.0307 |-0.0333 | 0.0440 |-0.0088 |-0.0076 |-0.0097 | 0.0241 | 0.0532 | 0.0213 |-2.8667 | 0.5622 |-3.2518
  Transform number: 3                  |-0.0204 | 0.0023 |-0.0292 | 0.0158 | 0.0264 |-0.0083 |-0.0071 |-0.0135 | 0.0055 | 2.3881 |-1.3371 | 0.1937
  Transform number: 4                  |-0.0186 | 0.0733 | 0.0095 |-0.0460 | 0.0168 | 0.0087 | 0.0651 | 0.0341 |-0.0072 |-0.4086 | 5.3377 |-6.0357



Data to excel output
--------------------
   0 (0)|  1 (10)        | 2(15)|   3(9)     |      4(12)       |   5(11)          |
pdb code| num ncs copies | year | resolution | data completeness| solvent fraction |

 6(13)         |   7(14)
data/param ASU | data/param NCS |


 8(1)             |     9(2)          |     10(5)      | 11(6)
R-work pdb header | R-free pdb header | R-work init    | init

  12(7)(16)               |        13(8)(16)          |
R-work without strict ncs | R-free without strict ncs |

    14(7)(16)             |    15(8)(16)
 R-work using strict ncs | R-free using strict ncs |


      16(7)(18)        |  17(8)(18)           |
 R-work with transform |R-free with transform |

    18(17)
 use geometry restraints |

 19(19)(16)   | 18(19)(16) |    20(19)(18)  |
time no ncs   | time ncs   | time_transform |

"""

class results_collection(object):

  def __init__(self):
    self.data_records = []
    self.data_records_long = []
    self.data_records_strict_ncs_dict = {}
    self.data_records_transform_dict = {}
    self.data_records_without_strict_ncs_dict = {}
    self.data_records_dict = {}
    self.transform_records_dict = {}
    self.cols_inp = [
      'pdb_code',
      'r_work_pdb_header',
      'r_free_pdb_header',
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
      'data_to_atoms_ratio_asu',
      'data_to_atoms_ratio_ncs',
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
      'data_to_atoms_ratio_asu',
      'data_to_atoms_ratio_ncs',
      'r_work_pdb_header',
      'r_free_pdb_header',
      'r_work_asu_init',
      'r_free_asu_init',
      'r_work_asu_no_ncs',
      'r_free_asu_no_ncs',
      'r_work_asu_ncs',
      'r_free_asu_ncs',
      'r_work_asu_transform',
      'r_free_asu_transform',
      'use_geometry_restraints',
      'time_no_ncs',
      'time_ncs',
      'time_transform']

    self.cols_names_for_table1 = [x.replace('_',' ') for x in self.cols_inp]
    self.cols_names_for_table2 = [x.replace('_',' ') for x in self.cols_names]
    # map the location in cols_inp to cols_names
    self.map_to_no_ncs = self.get_map(use_strict_ncs=False)
    self.map_to_ncs = self.get_map(use_strict_ncs=True,use_transforms=False)
    self.map_to_transform=self.get_map(use_strict_ncs=True,use_transforms=True)
    # set records type
    self.set_records_type()

    class sort_inp():
      def __init__(self):
        """ Parameter number to sort by """
    sort_inp.__dict__.update(zip(self.cols_inp, range(len(self.cols_inp))))
    self.sort_inp = sort_inp

    class sort_out():
      def __init__(self):
        """ Parameter number to sort by """
    sort_out.__dict__.update(zip(self.cols_names, range(len(self.cols_names))))
    self.sort_out = sort_out

  def set_records_type(self):
    """ Create list containing the type of each record """
    strings = ['pdb_code']
    floats = [
      'r_work_pdb_header',
      'r_free_pdb_header',
      'r_work_single_ncs',
      'r_free_single_ncs',
      'r_work_asu_init',
      'r_free_asu_init',
      'r_work_asu_final',
      'r_free_asu_final',
      'resolution',
      'solvent_fraction',
      'data_completeness',
      'data_to_atoms_ratio_asu',
      'data_to_atoms_ratio_ncs',
      'time']
    integers = ['num_ncs_copies','year']
    bools = ['use_strict_ncs', 'use_geometry_restraints', 'use_transforms']

    r = range(len(self.cols_inp))
    self.float_type_records = [i for i in r if self.cols_inp[i] in floats]
    self.int_type_records = [i for i in r if self.cols_inp[i] in integers]
    self.bool_type_records = [i for i in r if self.cols_inp[i] in bools]
    self.str_type_records = [i for i in r if self.cols_inp[i] in strings]

  def get_map(self,use_strict_ncs=True,use_transforms=False):
    """
    set the map between the input records and the combined output line

    Args:
      use_strict_ncs (bool): when True, using strict NCS
      use_transforms (bool): when True, refining rotation and translation

    Returns:
      (dict): map record location in the input to the record number at the
        output
    """
    if not use_strict_ncs:
      d = 0
    elif use_transforms:
      d = 2
    else:
      d = 1

    direct_match = ['pdb_code','r_work_pdb_header','r_free_pdb_header',
                    'r_work_asu_init','r_free_asu_init','resolution',
                    'num_ncs_copies','solvent_fraction','data_completeness',
                    'data_to_atoms_ratio_asu','data_to_atoms_ratio_ncs','year',
                    'use_geometry_restraints']
    type_match = [['time',['time_no_ncs','time_ncs','time_transform']],
                  ['r_work_asu_final',['r_work_asu_no_ncs','r_work_asu_ncs',
                                       'r_work_asu_transform',]],
                  ['r_free_asu_final',['r_free_asu_no_ncs','r_free_asu_ncs',
                                       'r_free_asu_transform',]]]
    map_dict = {}
    r1 = range(len(self.cols_inp))
    r2 = range(len(self.cols_names))
    for val in direct_match:
      k = [i for i in r1 if val == self.cols_inp[i]][0]
      v = [i for i in r2 if val == self.cols_names[i]][0]
      map_dict[k] = v

    for v1,v_all in type_match:
      v2 = v_all[d]
      k = [i for i in r1 if v1 == self.cols_inp[i]][0]
      v = [i for i in r2 if v2 == self.cols_names[i]][0]
      map_dict[k] = v

    return map_dict


  def read_filenames(self):
    """
    """
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
        print_warning(ln,pdb_code)
        # process files with good data (line start with pdb code)
        if re.match(regex,ln[:13]):
          data = self.process_data_according_to_data_type(ln)
          # collect data as is
          self.data_records.append(data)
          # combine all data for each pdb file to one row of data
          data_record = self.combine_records(data)
          self.data_records_dict[data[0]] = data_record
          # build other dictionaries
          if data[use_rec.use_strict_ncs]:
            self.data_records_strict_ncs_dict[data[0]] = data
          else:
            self.data_records_without_strict_ncs_dict[data[0]] = data
          # build other dictionaries
          if data[use_rec.use_transforms]:
            self.data_records_transform_dict[data[0]] = data
        # process transformation info
        if ln.startswith('  Transform number:'):
          transform_data = [x.strip() for x in ln.split('|')][1:]
          transform_data = map(float,transform_data)
          if self.transform_records_dict.has_key(data[0]):
            self.transform_records_dict[data[0]].append(transform_data)
          else:
            self.transform_records_dict[data[0]] = [transform_data]

    # create a list version of the dictionary
    for key,val in self.data_records_dict.iteritems():
      self.data_records_long.append(val)
    print '-'*60
    print 'There are {} good data records'.format(len(self.data_records))
    print 'There are {} pdb records'.format(len(self.data_records_long))
    print '-'*60

    # organize data
    # Sort records
    # self.data_records.sort(key=lambda x:x[use_rec.time])
    print 'done with collection'
    pickle.dump(self.data_records,open('data_records_v2','w'))
    pickle.dump(self.data_records_long,open('data_records_long_v2','w'))
    pickle.dump(self.data_records_dict,open('data_records_dict_v2','w'))
    pickle.dump(self.transform_records_dict,open('transform_records_dict','w'))
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
    self.transform_records_dict =pickle.load(open('transform_records_dict','r'))
    self.data_records_strict_ncs_dict = pickle.load(
      open('data_records_strict_ncs_dict_v2','r'))
    self.data_records_without_strict_ncs_dict = pickle.load(
      open('data_records_without_strict_ncs_dict_v2','r'))
    print 'Got the data from files'

  def plot_results_1(self):
    """
    Plot a result summery:
    r-work (refined using strict-ncs) vs.
    r-work (refined without strict-ncs)

    ncs_err , no_ncs_err are the difference between r-work and r-free for the two
    refinement methods

    The size of the circles indicate the size of (ncs_err - no_ncs_err). It is blue
    when the strict_ncs refinement make the (r_work - r_free) smaller.
    """
    rec_num = self.sort_inp
    # Collect data points that have results with and without strict_ncs
    x_ncs = []; y_no_ncs =[]
    x2 = []; y2 = []
    time_ncs = []; time_no_ncs = []
    pdb_code = []
    for x in self.data_records_strict_ncs_dict:
      if self.data_records_without_strict_ncs_dict.has_key(x):
        x_ncs.append(self.data_records_strict_ncs_dict[x]
                           [rec_num.r_work_asu_final])
        y_no_ncs.append(self.data_records_without_strict_ncs_dict[x]
                                  [rec_num.r_work_asu_final])
        x2.append(self.data_records_strict_ncs_dict[x]
                      [rec_num.r_free_asu_final])
        y2.append(self.data_records_without_strict_ncs_dict[x]
                      [rec_num.r_free_asu_final])
        time_ncs.append(self.data_records_strict_ncs_dict[x]
                      [rec_num.time])
        time_no_ncs.append(self.data_records_without_strict_ncs_dict[x]
                      [rec_num.time])
        pdb_code.append(x)

    # variable error bar values
    ncs_err = [abs(x-y) for (x,y) in zip(x_ncs,x2)]
    no_ncs_err = [abs(x-y) for (x,y) in zip(y_no_ncs,y2)]

    # get larges r-values for the 45 degrees line
    maxval = max(x_ncs + y_no_ncs) *1.05

    # Convert lists to numpy arrays
    x_ncs = np.array(x_ncs)
    y_no_ncs = np.array(y_no_ncs)
    ncs_err = np.array(ncs_err)
    no_ncs_err = np.array(no_ncs_err)

    # First illustrate basic pyplot interface, using defaults where possible.
    plt.figure()
    colors = []
    # delta_errs = [(x-y) for (x,y) in zip (x_ncs,y_no_ncs)]
    delta_errs = [(x-y) for (x,y) in zip (ncs_err,no_ncs_err)]
    # delta_errs_percent=[(x-y)/r for (x,y,r) in zip (ncs_err,no_ncs_err,x_ncs)]
    point_size = lambda x: 5 + 1000 * abs(x)
    s = [point_size(x) for x in delta_errs]
    for i,delta_err in enumerate(delta_errs):
      if delta_err > 0:
        colors.append('y')
        print 'Smaller (r_work - r_free) without NCS: ',pdb_code[i]
      else:
        colors.append('b')

    # add points for size reference : Difference of 0, 0.05, 0.1
    s.extend([point_size(0.1),point_size(0.05),point_size(0)])
    colors.extend(['g','g','g'])
    x_ncs = np.append(x_ncs,[.02,.02,.02])
    d = maxval - 0.55
    y_no_ncs = np.append(y_no_ncs,[.4 + d,.45 + d,.5 + d])
    plt.scatter(x_ncs,y_no_ncs,s,c=colors)
    plt.plot([0,maxval],[0,maxval])
    plt.xlabel('r-work Refinement using strict-ncs')
    plt.ylabel('r-work Refinement without strict-ncs')
    plt.title('strict-ncs influence on (R_work - R_free)')
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

    # print time outliers
    # d = 0.1
    # for (x,y,code) in zip(time_ncs,time_no_ncs,pdb_code):
    #   if abs(x-y)/x > d:
    #     print 'Time difference is {0:.0f}% for {1}'.format(100*(x-y)/x,code)

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

  def plot_results_3(self):
    """
    Plot a result summery:
    r-work (refined using strict-ncs and Transforms) vs.
    r-work (refined without strict-ncs)

    ncs_err , no_ncs_err are the difference between r-work and r-free for the two
    refinement methods

    The size of the circles indicate the size of (ncs_err - no_ncs_err). It is blue
    when the strict_ncs refinement make the (r_work - r_free) smaller.
    """
    sort_by = self.sort_inp
    # Collect data points that have results with and without strict_ncs
    x_ncs = []; y_no_ncs =[]
    x2 = []; y2 = []
    pdb_code = []
    print len(self.data_records_transform_dict)
    for x in self.data_records_transform_dict:
      if self.data_records_without_strict_ncs_dict.has_key(x):
        x_ncs.append(self.data_records_transform_dict[x]
                           [sort_by.r_work_asu_final])
        y_no_ncs.append(self.data_records_without_strict_ncs_dict[x]
                                  [sort_by.r_work_asu_final])
        x2.append(self.data_records_transform_dict[x]
                      [sort_by.r_free_asu_final])
        y2.append(self.data_records_without_strict_ncs_dict[x]
                      [sort_by.r_free_asu_final])
        pdb_code.append(x)

    # get larges r-values for the 45 degrees line
    maxval = max(x_ncs + y_no_ncs) *1.05

    # variable error bar values
    ncs_err = [abs(x-y) for (x,y) in zip(x_ncs,x2)]
    no_ncs_err = [abs(x-y) for (x,y) in zip(y_no_ncs,y2)]

    # Convert lists to numpy arrays
    x_ncs = np.array(x_ncs)
    y_no_ncs = np.array(y_no_ncs)

    # First illustrate basic pyplot interface, using defaults where possible.
    plt.figure()
    colors = []
    # delta_errs = [(x-y) for (x,y) in zip (x_ncs,y_no_ncs)]
    delta_errs = [(x-y) for (x,y) in zip (ncs_err,no_ncs_err)]

    point_size = lambda x: 5 + 1000 * abs(x)
    s = [point_size(x) for x in delta_errs]
    for i,delta_err in enumerate(delta_errs):
      if x > 0:
        colors.append('y')
        print 'Smaller (r_work - r_free) without NCS: ',pdb_code[i]
      else:
        colors.append('b')

    # add points for size reference : Difference of 0, 0.05, 0.1
    s.extend([point_size(0.1),point_size(0.05),point_size(0)])
    colors.extend(['g','g','g'])
    x_ncs = np.append(x_ncs,[.02,.02,.02])
    d = maxval - 0.55
    y_no_ncs = np.append(y_no_ncs,[.4 + d,.45 + d,.5 + d])
    plt.scatter(x_ncs,y_no_ncs,s,c=colors)
    plt.plot([0,maxval],[0,maxval])
    plt.xlabel('R_work Refinement, strict-ncs and transforms',fontsize=14)
    plt.ylabel('R_work Refinement, without strict-ncs',fontsize=14)
    plt.title('Transform and strict-NCS influence on (R_work - R_free)',fontsize=14)
    # plot reference point, to indicate the meaning of size
    plt.text(0.035,0.49 + d, 'The same r_work - r_free value',fontsize=14)
    plt.text(0.035,0.44 + d, '0.05 difference',fontsize=14)
    plt.text(0.035,0.39 + d, '0.1 difference',fontsize=14)

    plt.xlim(0,maxval)
    plt.ylim(0,maxval)
    plt.savefig('transform_refinement.png',transparent=False)
    plt.show()

  def plot_results_4(self):
    """
    Plot a result summery:
    r-free (refined using strict-ncs and Transforms) vs.
    r-free (refined without strict-ncs)
    """
    s = self.sort_inp
    # Collect data points that have results with and without strict_ncs
    x_ncs = []; y_no_ncs =[]
    pdb_code = []
    print len(self.data_records_transform_dict)
    for x in self.data_records_transform_dict:
      if self.data_records_without_strict_ncs_dict.has_key(x):
        x_ncs.append(
          self.data_records_transform_dict[x][s.r_free_asu_final])
        y_no_ncs.append(
          self.data_records_without_strict_ncs_dict[x][s.r_free_asu_final])
        pdb_code.append(x)

    # get larges r-values for the 45 degrees line
    maxval = max(x_ncs + y_no_ncs) *1.05

    # Convert lists to numpy arrays
    x_ncs = np.array(x_ncs)
    y_no_ncs = np.array(y_no_ncs)

    # First illustrate basic pyplot interface, using defaults where possible.
    plt.figure()

    plt.plot([0,maxval],[0,maxval])
    plt.plot(x_ncs,y_no_ncs,'o')
    plt.xlabel('With strict-ncs and transforms',fontsize=14)
    plt.ylabel('Without strict-ncs',fontsize=14)
    plt.title('R-Free after refinement',fontsize=14)

    plt.xlim(0,maxval)
    plt.ylim(0,maxval)
    plt.savefig('r_free.png',transparent=False)
    plt.show()

  def plot_results_5(self):
    """
    Plot a result summery:
    r_free - r_work (refined using strict-ncs and Transforms) vs.
    r_free - r_work (refined without strict-ncs)
    """
    s = self.sort_inp
    # Collect data points that have results with and without strict_ncs
    x_ncs = []; y_no_ncs =[]
    x2 = []; y2 = []
    pdb_code = []
    print len(self.data_records_transform_dict)
    for x in self.data_records_transform_dict:
      if self.data_records_without_strict_ncs_dict.has_key(x):
        x_ncs.append(self.data_records_strict_ncs_dict[x]
                           [s.r_work_asu_final])
        y_no_ncs.append(self.data_records_without_strict_ncs_dict[x]
                                  [s.r_work_asu_final])
        x2.append(self.data_records_strict_ncs_dict[x]
                      [s.r_free_asu_final])
        y2.append(self.data_records_without_strict_ncs_dict[x]
                      [s.r_free_asu_final])
        pdb_code.append(x)

    # variable error bar values
    ncs_err = [abs(x-y) for (x,y) in zip(x_ncs,x2)]
    no_ncs_err = [abs(x-y) for (x,y) in zip(y_no_ncs,y2)]

    # Convert lists to numpy arrays
    x_ncs = np.array(ncs_err)
    y_no_ncs = np.array(no_ncs_err)

    # get larges r-values for the 45 degrees line
    maxval = max(x_ncs + y_no_ncs) *1.05
    maxval = 0.2

    # First illustrate basic pyplot interface, using defaults where possible.
    plt.figure()

    plt.plot([0,maxval],[0,maxval])
    plt.plot(x_ncs,y_no_ncs,'o')
    plt.xlabel('With strict-ncs and transforms',fontsize=14)
    plt.ylabel('Without strict-ncs',fontsize=14)
    plt.title('(R_Free - R_work) after refinement',fontsize=14)

    plt.xlim(0,maxval)
    plt.ylim(0,maxval)
    plt.savefig('r_free-r_work.png',transparent=False)
    plt.show()

  def plot_results_6(self):
    """
    Plot a result summery:
    r-free (refined using strict-ncs) vs.
    r-free (refined without strict-ncs)

    ncs_err , no_ncs_err are the difference between r-work and r-free for the two
    refinement methods

    The size of the circles indicate the size of (ncs_err - no_ncs_err). It is blue
    when the strict_ncs refinement make the (r_work - r_free) smaller.
    """
    rec_num = self.sort_inp
    # Collect data points that have results with and without strict_ncs
    x_ncs = []; y_no_ncs =[]
    x2 = []; y2 = []
    time_ncs = []; time_no_ncs = []
    pdb_code = []
    for x in self.data_records_strict_ncs_dict:
      if self.data_records_without_strict_ncs_dict.has_key(x):
        x_ncs.append(self.data_records_strict_ncs_dict[x]
                           [rec_num.r_work_asu_final])
        y_no_ncs.append(self.data_records_without_strict_ncs_dict[x]
                                  [rec_num.r_work_asu_final])
        x2.append(self.data_records_strict_ncs_dict[x]
                      [rec_num.r_free_asu_final])
        y2.append(self.data_records_without_strict_ncs_dict[x]
                      [rec_num.r_free_asu_final])
        time_ncs.append(self.data_records_strict_ncs_dict[x]
                      [rec_num.time])
        time_no_ncs.append(self.data_records_without_strict_ncs_dict[x]
                      [rec_num.time])
        pdb_code.append(x)

    # variable error bar values
    ncs_err = [abs(x-y) for (x,y) in zip(x_ncs,x2)]
    no_ncs_err = [abs(x-y) for (x,y) in zip(y_no_ncs,y2)]

    # get larges r-values for the 45 degrees line
    maxval = max(x2 + y2) *1.05

    # Convert lists to numpy arrays
    x2 = np.array(x2)
    y2 = np.array(y2)
    ncs_err = np.array(ncs_err)
    no_ncs_err = np.array(no_ncs_err)

    # First illustrate basic pyplot interface, using defaults where possible.
    plt.figure()
    colors = []
    # delta_errs = [(x-y) for (x,y) in zip (x_ncs,y_no_ncs)]
    delta_errs = [(x-y) for (x,y) in zip (ncs_err,no_ncs_err)]
    point_size = lambda x: 5 + 1000 * abs(x)
    s = [point_size(x) for x in delta_errs]
    for i,delta_err in enumerate(delta_errs):
      if delta_err > 0:
        colors.append('y')
        print 'Smaller (r_work - r_free) without NCS: ',pdb_code[i]
      else:
        colors.append('b')

    # add points for size reference : Difference of 0, 0.05, 0.1
    s.extend([point_size(0.1),point_size(0.05),point_size(0)])
    colors.extend(['g','g','g'])
    x2 = np.append(x2,[.02,.02,.02])
    d = maxval - 0.55
    p1_y = 0.48
    dy = 0.03
    dr = 0.01
    y2 = np.append(y2,[p1_y - 2*dy + d,p1_y - dy + d,p1_y + d])
    plt.scatter(x2,y2,s,c=colors)
    plt.plot([0,maxval],[0,maxval])
    plt.xlabel('R-free with strict-ncs',fontsize=14)
    plt.ylabel('R-free without strict-ncs',fontsize=14)
    plt.title('strict-ncs influence on (R_work - R_free)',fontsize=14)
    # plot reference point, to indicate the meaning of size
    plt.text(0.010,p1_y + dy + d, '(R_free - R_work) difference',fontsize=16)
    plt.text(0.035,p1_y + d - dr, '0.00',fontsize=14)
    plt.text(0.035,p1_y - dy + d - dr, '0.05 ',fontsize=14)
    plt.text(0.035,p1_y - 2*dy + d - dr, '0.10 ',fontsize=14)

    plt.xlim(0,maxval)
    plt.ylim(0,maxval)
    plt.savefig('r_free_bubble_plot.png',transparent=False)
    plt.show()

  def save_csv_table_to_file(self):
    """
    Save a the collected data in a csv file using the record order as
    specified in self.cols_names_for_table1
    """
    print 'Running ',sys._getframe().f_code.co_name
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

    The data in the file is one line per test with the records order as
    specified in self.cols_names_for_table2
    """
    print 'Running ',sys._getframe().f_code.co_name
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

  def save_csv_selected_files(self,pdb_file_list):
    """
    Create a CSV of selected files in a modified format
    """
    print 'Running ',sys._getframe().f_code.co_name
    ur = self.sort_out
    file_name = 'selected_files_v2.csv'
    records_to_use = self.cols_names_for_table2[:ur.r_free_asu_transform]
    table_title = ['PDB','n','Year','Res.','Comp','Solvent',
                   'D/P ASU','D/P NCS',
                   'Header','Initial','W/O NCS','W NCS','Transform']

    f = open(file_name,'w')
    f.write(','.join(table_title))
    f.write('\n')
    for rec in self.data_records_long:
      if rec[ur.pdb_code] in pdb_file_list:
        rec = [str(x) for x in rec]
        outstr = rec[0:8]
        outstr += ['/'.join(rec[ur.r_work_pdb_header:ur.r_free_pdb_header+1]),
                  '/'.join(rec[ur.r_work_asu_init:ur.r_free_asu_init+1]),
                  '/'.join(rec[ur.r_work_asu_no_ncs:ur.r_free_asu_no_ncs+1]),
                  '/'.join(rec[ur.r_work_asu_ncs:ur.r_free_asu_ncs+1]),
                  '/'.join(rec[ur.r_work_asu_transform:ur.r_free_asu_transform+1])]
        outstr += '\n'
        f.write(','.join(outstr))
    f.close()
    print 'data was saved to: ',file_name

  def get_list_of_unprocessed_files(self):
    """
    Check which of the PDB files, from our initial list, are not included in
    the results
    """
    print 'Running ',sys._getframe().f_code.co_name
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
    with_transform_set=set(self.data_records_transform_dict.keys())
    print '\n  files with problems'
    print '-----------------------------------'
    outdata1 = pdb_code_set - without_strict_ncs_set
    outdata2 = pdb_code_set - with_strict_ncs_set
    outdata3 = pdb_code_set - with_transform_set
    outdata = outdata3.union(outdata2)
    outdata = list(outdata.union(outdata1))

    n = 5
    l = len(outdata)
    print 'Total number of files: ',len(pdb_code_set)
    print(list(outdata))
    for i in range(0,l,n):
      e = min(i+n,l)
      s = ['{}']*(e-i)
      s = ', '.join(s)
      print s.format(*outdata[i:e])

    print '\nGood files: ', len(pdb_code_set) - len(outdata)

  def filter_best_1(self):
    """
    Filter files with best (r_work_asu_init - r_work_asu_transform)
    """
    print '\nFilter files with best (r_work_asu_init - r_work_asu_transform)'
    print '---------------------------------------------------------------'
    temp = []
    file_list = []
    use_rec = self.sort_out
    for rec in self.data_records_long:
      x2,x1,use_data = self.final_vs_init(rec)
      if use_data:
        if rec[use_rec.r_work_asu_transform] != \
                rec[use_rec.r_free_asu_transform]:
          temp.extend([[round(x1-x2,4),rec]])
    temp.sort()
    for rec in temp[-15:]:
      print rec
      file_list.append(rec[1][0])
    print set([x[1][0] for x in temp])

  def filter_best_2(self):
    """
    Filter files with r_free_asu_transform == r_work_asu_transform
    """
    print '\nFilter files with r_free_asu_transform == r_work_asu_transform'
    print '--------------------------------------------------------------'
    temp = []
    file_list = []
    use_rec = self.sort_out
    for rec in self.data_records_long:
      x1 =  rec[use_rec.r_free_asu_transform]
      y1 =  rec[use_rec.r_work_asu_transform]
      if  isinstance(x1,float):
        if x1 == y1:
          temp.extend([[x1,rec]])
    temp.sort()
    for rec in temp[-10:]:
      file_list.append(rec[1][0])
      print rec
    print set([x[1][0] for x in temp])

  def filter_best_3(self):
    """
    Filter files with best delta(r_work - r_free) improvement
    """
    print '\nFilter files with best delta(r_work - r_free) improvement'
    print '---------------------------------------------------------'
    temp = []
    file_list = []
    use_rec = self.sort_out
    for rec in self.data_records_long:
      x1 =  rec[use_rec.r_work_asu_ncs]
      x2 =  rec[use_rec.r_free_asu_ncs]
      y1 =  rec[use_rec.r_work_asu_transform]
      y2 =  rec[use_rec.r_free_asu_transform]
      if  isinstance(x1,float) and isinstance(y1,float):
        temp.extend([[round(abs(x1-x2) - abs(y1-y2),4),rec]])
    temp.sort()
    for rec in temp[-10:]:
      file_list.append(rec[1][0])
      print rec
    print set([x[1][0] for x in temp])

  def filter_best_4(self):
    """
    Collect files with final R-work larger than initial
    """
    print '\nFilter files with final R-work larger than initial'
    print '--------------------------------------------------'
    temp = []
    file_list = []
    use_rec = self.sort_out
    for rec in self.data_records_long:
      x1 =  rec[use_rec.r_work_asu_init]
      y1 =  rec[use_rec.r_work_asu_transform]
      if  isinstance(x1,float) and isinstance(y1,float):
        if y1 > x1:
          temp.extend([[round(y1-x1,4),rec]])
    temp.sort()
    for rec in temp[-10:]:
      file_list.append(rec[1][0])
      print rec
    print set([x[1][0] for x in temp])
    print 'number of files: ',len(temp)

  def r_free_better_without_ncs(self):
    """
    Collect/show structures where the r_free_better_without_ncs

    Argument:
    rec : list of all refinement data
    Return:
    f,i : (float,float) R-free: final,initial
    use_data : (bool) When True, both values are floats
    """
    print '\nr_free_better_without_ncs'
    print '--------------------------------------------------'
    temp = []
    file_list = []
    use_rec = self.sort_out
    for rec in self.data_records_long:
      x1 =  rec[use_rec.r_free_asu_no_ncs]
      y1 =  rec[use_rec.r_free_asu_transform]
      if  isinstance(x1,float) and isinstance(y1,float):
        if y1 > x1:
          temp.extend([[round(y1-x1,4),rec]])
    temp.sort()
    for rec in temp[-10:]:
      file_list.append(rec[1][0])
      print rec
    print set([x[1][0] for x in temp])
    print 'number of files: ',len(temp)

  def big_change_in_transforms(self,print_all=False, print_n=5):
    """
    Collect files with big change in transformations
    """
    print '\nCollect files with big change in transformations'
    print '------------------------------------------------'
    list_of_files = []
    set_of_files = set()
    for k,v in self.transform_records_dict.iteritems():
      for rec in v:
        m = max(rec)
        if m > 0.2:
          list_of_files.append([m,self.data_records_dict[k]])
          set_of_files.add(k)
    list_of_files.sort()
    list_of_files.reverse()
    printed_files = []
    print 'number of files with large change in transforms: ',len(list_of_files)
    if print_all: print_n = len(list_of_files)
    i = 0
    for l in list_of_files:
      if i == print_n: break
      if l[1][0] not in printed_files:
        print l
        i += 1
        printed_files.append(l[1][0])
    print printed_files[:print_n]

  def filter_best_6(self):
    """
    Collect files with R_work (initial - Header) > 0.2
    """
    print '\nCollect files with R_work (initial - Header) > 0.2'
    print '--------------------------------------------------'
    temp = []
    file_list = []
    use_rec = self.sort_out
    for rec in self.data_records_long:
      x1 =  rec[use_rec.r_work_asu_init]
      y1 =  rec[use_rec.r_work_pdb_header]
      if  isinstance(x1,float) and isinstance(y1,float):
        if (y1 != -1) and (x1 - y1) > 0.2:
          temp.extend([[round(abs(y1-x1),4),rec]])
    temp.sort()
    for rec in temp:
      file_list.append(rec[1][0])
      print rec
    print set([x[1][0] for x in temp])
    print 'number of files: ',len(temp)

  def filter_best_7(self):
    """
    Collect files with R_work final > 0.5
    """
    print '\nCollect files with R_work final > 0.5'
    print '-------------------------------------'
    temp = []
    file_list = []
    use_rec = self.sort_out
    for rec in self.data_records_long:
      x1 =  rec[use_rec.r_work_asu_transform]
      if  isinstance(x1,float):
        if x1 > 0.5:
          temp.extend([[round(x1,4),rec]])
    temp.sort()
    for rec in temp:
      file_list.append(rec[1][0])
      print rec
    print set([x[1][0] for x in temp])
    print 'number of files: ',len(temp)

  def get_stat(self,exclude_list=[]):
    """
    Get mean and standard deviation of:
    1) improvement in R-work, All refinements methods compare to initial value
    2) improvement in delta(R_work - R_free)  All refinements methods compare to without NCS
    3) improvement in R-work, All refinements methods compare to without transform

    exclude_list: (list of str) contains list of outliers not to include
    """
    use_rec = self.sort_out
    st1 ='R-work improvement, All refinements methods compare to initial value'
    sf1 = self.final_vs_init
    st2 ='delta(R_work - R_free) improvement, All refinements methods compare to without NCS'
    sf2 = self.delta_work_free
    st3 ='R-work improvement, applying transform refinement'
    sf3 = self.final_vs_ncs
    st4 = 'delta(R_work - R_free) All refinements methods compare to PDB'
    sf4 = self.delta_work_free_pdb
    get_stat_using = [(st1,sf1),(st2,sf2),(st3,sf3),(st4,sf4)]
    for test in get_stat_using:
      i = []
      f = []
      print test[0]
      for rec in self.data_records_long:
        f1,i1,use_data =test[1](rec)
        if use_data:
          i.append(i1)
          f.append(f1)
      i = np.array(i)
      f = np.array(f)
      d = i - f
      print 'Mean : {0:.4f}  STD : {1:.4f}   Number of files {2}'\
        .format(np.mean(d),np.std(d),len(d))

  def final_vs_init(self,rec):
    """
    Compare R-work, All refinements methods compare to initial value

    Argument:
    rec : list of all refinement data
    Return:
    f,i : (float,float) R-work : final,initial
    use_data : (bool) When True, both values are floats
    """
    use_rec = self.sort_out
    f =  rec[use_rec.r_work_asu_transform]
    i =  rec[use_rec.r_work_asu_init]
    use_data = isinstance(f,float) and isinstance(i,float)
    return f,i,use_data

  def final_vs_pdb(self,rec):
    """
    Compare R-work, All refinements methods compare to pdb published value

    Argument:
    rec : list of all refinement data
    Return:
    f,i : (float,float) R-work : final,initial
    use_data : (bool) When True, both values are floats
    """
    use_rec = self.sort_out
    f =  rec[use_rec.r_work_asu_transform]
    i =  rec[use_rec.r_work_pdb_header]
    use_data = isinstance(f,float) and isinstance(i,float)
    return f,i,use_data

  def delta_work_free(self,rec):
    """
    Compare delta(R_work - R_free)  All refinements methods compare to without NCS

    Argument:
    rec : list of all refinement data
    Return:
    f,i : (float,float) delta(R_work - R_free) : final,initial
    use_data : (bool) When True, both values are floats
    """
    use_rec = self.sort_out
    i1 =  rec[use_rec.r_work_asu_no_ncs]
    i2 =  rec[use_rec.r_free_asu_no_ncs]
    f1 =  rec[use_rec.r_work_asu_transform]
    f2 =  rec[use_rec.r_free_asu_transform]
    use_data = isinstance(f1,float) and isinstance(i1,float)
    if use_data:
      f =  abs(f2-f1)
      i =  abs(i2-i1)
    else:
      f = i = None
    return f,i,use_data

  def delta_work_free_pdb(self,rec):
    """
    Compare delta(R_work - R_free)  All refinements methods compare to PDB

    Argument:
    rec : list of all refinement data
    Return:
    f,i : (float,float) delta(R_work - R_free) : final,initial
    use_data : (bool) When True, both values are floats
    """
    use_rec = self.sort_out
    i1 =  rec[use_rec.r_work_pdb_header]
    i2 =  rec[use_rec.r_free_pdb_header]
    f1 =  rec[use_rec.r_work_asu_transform]
    f2 =  rec[use_rec.r_free_asu_transform]
    test1 = isinstance(i1,float) and isinstance(i2,float)
    test2 = test1 and (-1 not in [i1,i2])
    use_data = isinstance(f1,float) and test2
    if use_data:
      f =  abs(f2-f1)
      i =  abs(i2-i1)
    else:
      f = i = None
    return f,i,use_data

  def final_vs_ncs(self,rec):
    """
    Compare R-work, All refinements methods compare to without transform

    Argument:
    rec : list of all refinement data
    Return:
    f,i : (float,float) R-work : final,initial
    use_data : (bool) When True, both values are floats
    """
    use_rec = self.sort_out
    f =  rec[use_rec.r_work_asu_transform]
    i =  rec[use_rec.r_work_asu_ncs]
    use_data = isinstance(f,float) and isinstance(i,float)
    return f,i,use_data

  def filter_by_year(self,rec_list,year):
    """d
    filter a list of [[float,list],[float,list]...] by the year value in list
    Return rec_list with records from 'year' and on
    """
    use_rec = self.sort_out
    return [x for x in rec_list if x[1][use_rec.year] >= year]

  def process_data_according_to_data_type(self,ln):
    """ Process line of data from data files  """
    data = [x.strip() for x in ln.split('|')]
    # convert data values to proper type
    for indx in self.float_type_records:
      data[indx] = float(data[indx])
    for indx in self.int_type_records:
      if data[indx] == 'None': data[indx] = 0
      else: data[indx] = int(data[indx])
    for indx in self.bool_type_records:
      data[indx] = (data[indx] == 'True')
    return data

  def combine_records(self,data):
    """
    Every pdb record can be refined with different options. Collect and
    combine all records related to a pdb record to a single data line
    """
    use_rec = self.sort_inp
    if self.data_records_dict.has_key(data[0]):
      data_record = self.data_records_dict[data[0]]
    else:
      data_record = len(self.cols_names)*[None,]
    if data[use_rec.use_strict_ncs] and not data[use_rec.use_transforms]:
      for map_from,map_to in self.map_to_ncs.iteritems():
        if not  data_record[map_to]:
          data_record[map_to] = data[map_from]
    if not data[use_rec.use_strict_ncs] and not data[use_rec.use_transforms]:
      for map_from,map_to in self.map_to_no_ncs.iteritems():
        if not  data_record[map_to]:
          data_record[map_to] = data[map_from]
    if data[use_rec.use_strict_ncs] and data[use_rec.use_transforms]:
      for map_from,map_to in self.map_to_transform.iteritems():
        if not  data_record[map_to]:
          data_record[map_to] = data[map_from]
    return data_record

def print_warning(ln,pdb_code):
  if ln.startswith('Warning'):
    print '{0}: {1}'.format(pdb_code,ln.strip())
  if ln.startswith('Using pdb file from local machine'):
    print 'Used {} from local machine'.format(pdb_code)
  if ln.startswith('Traceback'):
    print '{0}: {1}'.format(pdb_code,'Problem processing')


if __name__=='__main__':
  # fn = r'C:\Phenix\Dev\Work\work\NCS\junk\pdb_test\Refinement_data_9_12_2014'
  # fn = "/net/cci/youval/Work/work/NCS/junk/pdb_test/queue_job"
  fn = r"C:\Phenix\Dev\Work\work\NCS\junk\pdb_test\queue_job"

  path_to_log_files = fn
  current_path = os.getcwd()
  os.chdir(path_to_log_files)
  process_results = results_collection()

  # need to read_filenames only for new runs, it overwrites existing data files
  process_results.read_filenames()

  # process existing files
  process_results.get_data_from_files()
  process_results.save_csv_table_to_file()
  process_results.save_csv_table2_to_file()
  print '-'*30

  # plot
  process_results.plot_results_1()
  process_results.plot_results_2()
  process_results.plot_results_3()
  process_results.plot_results_4()
  process_results.plot_results_5()
  process_results.plot_results_6()

  # check which files are not processed yet
  process_results.get_list_of_unprocessed_files()

  s1 = ['3raa','2bfu','3qpr','3ntt','1dzl','3s4g','1ei7',
        '2qij','2ztn','3nou','3hag','2g34']
  s2 = ['3raa','2bfu','3qpr','3ntt','1dzl','3s4g','1ei7',
        '2qij','3nou','3hag','2g34',
        '1dwn','1qjy', '1qju', '1a37', '2buk','1vcr','1tnv',
        '3r0r','1vb2','4gbt','2wws','2wff']


  process_results.filter_best_1()
  process_results.filter_best_2()
  process_results.filter_best_3()
  process_results.filter_best_4()
  process_results.big_change_in_transforms()
  process_results.r_free_better_without_ncs()
  process_results.filter_best_6()
  process_results.filter_best_7()
  process_results.save_csv_selected_files(s2)
  process_results.get_stat(exclude_list=[])

  os.chdir(current_path)
  print 'Done...'
