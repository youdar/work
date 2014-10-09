from __future__ import division
import cPickle as pickle
import os

def read_raw_data(save_to_file=False,file_name=''):
  """ collect data from log files  """
  print 'Start collecting data'
  print '-'*30
  # Read files in directory_path
  path = '/net/cci-filer2/raid1/home/youval/Work/work/NCS/pdb_surveys'
  data_path = path + '/pdb_survey_logs'
  files = os.listdir(data_path)
  files = [x for x in files if x.startswith('log_')]
  print 'Number of log files: ',len(files)
  data = []
  improper_rotations = []
  master_ncs_relation_issues = []
  for fn in files:
    d = open(os.path.join(data_path, fn), "r").readlines()
    if len(d) == 1:
      rd = d[0].split(',')
      if len(rd) == 11:
        # format data
        # pdb_code, n_mtrix_rec, t_min, t_max, t_simple,
        d1 = [rd[0],int(rd[1]),int(rd[2]),int(rd[3]),int(rd[4])]
        # t_min_gruops, t_max_gruops, t_simple_gruops,
        d2 = [int(rd[5]),int(rd[6]),int(rd[7])]
        # t_min_time, t_max_time, t_simple_time
        d3 = [float(rd[8]),float(rd[9]),float(rd[10])]
        data.append(d1 + d2 + d3)
      else:
        if 'Sorry: Rotation matrices are not proper!' in d[0]:
          improper_rotations.append(fn[-4:])
        elif 'Master NCS and Copy are very poorly related' in d[0]:
          master_ncs_relation_issues.append(fn[-4:])
        else:
          print rfn[-4:], d[0]
  print 'Number of records collected: ',len(data)
  if save_to_file:
    fn = os.path.join(path,file_name)
    pickle.dump(data, open(fn,'w'))
    pickle.dump(improper_rotations, open(os.path.join(path,'improper_rotations'),'w'))
    pickle.dump(master_ncs_relation_issues, open(os.path.join(path,'master_ncs_relation_issues'),'w'))
    print 'Pickled data was saved to: "{}"'.format(file_name)
  return data

def run():
  file_name = 'ncs_in_pdb_survey_data'
  if os.path.isfile(file_name):
    print 'Using existing data file'
    data = pickle.load(open(file_name,'r'))
    improper_rotations = pickle.load(open('improper_rotations','r'))
    master_ncs_relation_issues = pickle.load(open('master_ncs_relation_issues','r'))
  else:
    print 'processing new data'
    data = read_raw_data(save_to_file=True,file_name=file_name)
  # format data
  # d[0]: pdb_code, d[1]: n_mtrix_rec, d[2]: t_min, d[3]: t_max, d[4]: t_simple,
  # d[5]: t_min_gruops, d[6]: t_max_gruops, d[7]: t_simple_gruops,
  # d[8]: t_min_time, d[9]: t_max_time, d[10]: t_simple_time
  print '='*30
  print 'processing records'
  print '-'*30
  print 'Number of records: ',len(data)
  print '\n Records order:'
  print '0: pdb_code, 1: n_mtrix_rec, 2: min, 3: max, 4: simple,'
  print '5: min_gruops, 6: max_gruops, 7: simple_gruops,'
  print '8: min_time,9: max_time, 10: simple_time\n'


  bad_rec = [d for d in data if d[4]==-1]
  print '\nNumber of records with "simple" processing issues: ',len(bad_rec)
  print_rec(bad_rec)

  bad_rec = [d for d in data if -1 in d[2:4]]
  print '\nNumber of records with "cctbx" processing issues: ',len(bad_rec)
  print_rec(bad_rec)

  d0 = [d for d in data if (d[1] == d[3])]
  print '\nNumber of records with equal max and MTRIX ncs operators: ',len(d0)
  print_rec(d0)

  d1 = [d for d in data if (d[1] != d[3])]
  print 'Number of records with different max and MTRIX ncs operators: ',len(d1)
  print_rec(d1)

  d1_1 = [d for d in data if ((d[1] ==0) and (d[3]==0))]
  s= 'Records with MTRIX ncs operators where ncs relation are present: '
  print s,len(d1_1)
  print_rec(d1_1)

  d2 = [d for d in data if (d[1] != d[2])]
  print 'Number of records with different min and MTRIX ncs operators: ',len(d2)
  print_rec(d2)

  d3 = [d for d in data if (d[3] != d[2])]
  print 'Number of records with different max and min ncs operators:   ',len(d3)
  print_rec(d3)

  d4 = [d for d in data if (d[3] != d[4])]
  print 'Number of records with different simple and min ncs operators:',len(d4)
  print_rec(d4)

  d5 = [d for d in data if (d[2] > d[1])]
  print 'Number of new ncs relation:                                   ',len(d5)
  print_rec(d5)

  d6 = [d for d in data if (d[4] in d[2:4])]
  print '"cctbx" find the number of ncs relations as "simple":         ',len(d6)
  print_rec(d6)

  d7 = [d for d in data if (d[0] in ['2wff','2wws'])]
  print '"2wws, 2wff":         ',len(d6)
  print_rec(d7)

  d_time = sorted(data, key= lambda x: x[9])
  print '10 slowest files to process: ', d_time[-10:]
  print_rec(d_time,reverse=True,nim_n=10,print_item=9)


  d_time = sorted(data, key= lambda x: x[10])
  print '10 slowest files to process: ', d_time[-10:]
  print_rec(d_time,reverse=True,nim_n=10,print_item=10)

  print 'Done...'

def print_rec(d,reverse=False,nim_n=5,print_item=None):
  n = min(nim_n,len(d))
  x = 1
  di = 0
  if reverse:
    x = -1
    di = 1
  for i in range(n):
    j = i * x - di
    if print_item:
      print d[j][print_item],d[j]
    else:
      print d[j]
  print '......'

if __name__ == '__main__':
  run()
