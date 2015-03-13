from __future__ import division
from random import shuffle
import collect_ncs_files
import cPickle as pickle
from glob import glob
import os


def run():
  c = collect_ncs_files.ncs_paper_data_collection()
  d = c.collect_all_file_records()
  print 'Total number of files in data dir:',len(d)
  file_list = glob(os.path.join(c.data_dir,'log_*'))
  file_list = [os.path.split(x)[-1] for x in file_list]
  file_list = {x[-4:] + '.pdb' for x in file_list}
  print len(file_list)
  ncs_only = {x for x in d if d[x].only_master_in_pdb}
  print 'number of files containing NCS only:',len(ncs_only)
  solvent_fraction = {x for x in d if d[x].solvent_fraction > 0.8}
  print 'number of files where solvent_fraction > 0.8:',len(solvent_fraction)
  completeness = {x for x in d if d[x].data_completeness > 0.8}
  print 'number of files where data_completeness > 0.8:',len(completeness)
  not_xray_diff = {x for x in d
                   if not ('X-RAY DIFFRACTION' in d[x].experiment_type)}
  print 'files not from X-RAY DIFFRACTION:',len(not_xray_diff)

  print '='*50
  asu_files = glob(os.path.join(c.asu_dir,'*.pdb'))
  asu_files = [os.path.split(x)[-1] for x in asu_files]
  asu_files = set(asu_files)
  print len(asu_files)
  print len(asu_files - file_list)
  print len(file_list - asu_files)

  print '='*50
  s = '/net/cci-filer2/raid1/home/youval/work/work/NCS/ncs_paper'
  s += '/ncs_paper_data_files/pdb_with_ncs_not_used'
  not_used_list = glob(os.path.join(s,'log_*'))
  not_used_list = [os.path.split(x)[-1] for x in not_used_list]
  not_used_list = {x[-4:] + '.pdb' for x in not_used_list}
  print len(not_used_list)
  print len(not_used_list - file_list)
  print len(file_list - not_used_list)


  print '='*50
  print "model_vs_data_files"
  model_vs_data_files = glob(os.path.join(c.model_vs_data_dir,'*.txt'))
  model_vs_data_files = [os.path.split(x)[-1] for x in model_vs_data_files]
  model_vs_data_files = {x[:4] + '.pdb' for x in model_vs_data_files}
  print len(model_vs_data_files)
  print len(asu_files - model_vs_data_files)
  print len(model_vs_data_files - asu_files)
  print asu_files - model_vs_data_files
  print model_vs_data_files - asu_files

  # print sample
  # keys = d.keys()
  # shuffle(keys)
  # keys = keys[:3]
  # keys = not_xray_diff
  # for k in keys:
  #   print '-'*50
  #   print d[k]


def look_at_files_not_used():
  collect = collect_ncs_files.ncs_paper_data_collection()
  d = {}
  files_list = glob(collect.pdb_not_used_dir + '/log_*')
  for fn in files_list:
    r = pickle.load(open(fn,'r'))
    d[r.pdb_id] = r
  #
  not_xray_diff = {x for x in d
                   if ('X-RAY DIFFRACTION' in d[x].experiment_type)}
  print 'files not from X-RAY DIFFRACTION:',len(not_xray_diff)


if __name__=='__main__':
  run()
  print '*'*50
  # look_at_files_not_used()
  print '*'*50
