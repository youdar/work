from __future__ import division
import cPickle as pickle
import os

'''
Collecting pdb files that have very different clashscore and nb_clashscore
So we can look at the cause and evaluate the difference in scoring methods
'''


  
def get_data():
  ''' () -> dict,dict,list
  Read results of clash score servey for PROBE clashscore and restraints manager nonbonded clashscore
  
  c:\Phenix\Dev\Work\work\Clashes\Data\clashscore_compare_ready_set-12-5-13_dict
  
  pdb_clash_scores = list([score_with_hydrogen,score_without_hydrogen]...)
  pdb_clash_score_and_name = list([score_with_hydrogen,score_without_hydrogen,experment_type,file_name]...)
  pdb_clash_score_dict[file_name] = [score_with_hydrogen,score_without_hydrogen,experment_type]
  
  Returns:
  data_dict: data_dict[pdb_file_name] = [total_nb_clashscore,without_sym_nb_clashscore,clashscore_probe]
  data_dict2: data_dict2[file_name] = [clashscore_all_clashes,clashscore_simple,clashscore_only_sym_op,clashscore_solvent_solvent]
  experiment_type_dict: experiment_type_dict[experiment_type] = list of pdb_file_name
  pdb_file_list: a list of all files that were compared (those are all the pdb files with clashscore < 50)
    
  
  >>>experiment_type_dict['NMR'][:10]
  ['103d', '124d', '141d', '142d', '169d', '175d', '1a1d', '1ac3', '1al5', '1anp']
  >>>data_dict['142d']
  [0.0, 0.0, 0.0]
  '''
  datapath = os.path.realpath('c:\Phenix\Dev\Work\work\Clashes\Data')
  data_dict_file = 'clashscore_compare_reduce_12_6_2013_dict' 	# Probe O vdw is 1.4
  #data_dict_file = 'clashscore_compare_reduce_12_11_2013_dict' 	# Probe O vdw is 1.52
  data_dict_file2 = 'clashscore_solvent_reduce_12_13_2013_dict' # with solvent-solvent results
  experiment_dict_file = 'experiment_type_to_files_dict'		# source for experiment_type_dict
  # Get data
  data_dict = pickle.load(open(os.path.join(datapath,data_dict_file),'r'))
  data_dict2 = pickle.load(open(os.path.join(datapath,data_dict_file2),'r'))
  experiment_type_dict = pickle.load(open(os.path.join(datapath,experiment_dict_file),'r'))
  # Collect all files that we compared 
  pdb_file_list = [key for key in data_dict]
  return data_dict,data_dict2,experiment_type_dict,pdb_file_list

def run():
  data_dict,data_dict2,experiment_type_dict,pdb_file_list = get_data()
  # Collect files where the abs(clashscore - nb_clashscore) > 10
  results = {}
  delta = 30
  delta_solvent = 10
  ratio = 0.3
  # Test data integrity: that old total score equals the new one
  for experiment_type in experiment_type_dict:
    for x in experiment_type_dict[experiment_type]:
      if data_dict2.has_key(x) and data_dict.has_key(x):
        if data_dict2[x][0] != data_dict[x][0]:
          print 'Clashscor totals do not match: ',x     
  #
  print 'The number of files with abs(clashscore - without_sym_nb_clashscore) > {}: '.format(delta)
  print '-'*85
  for experiment_type in experiment_type_dict:
    temp = [x for x in experiment_type_dict[experiment_type] 
            if data_dict2.has_key(x) and data_dict.has_key(x) and 
            (abs(data_dict2[x][1]-data_dict[x][2]) > delta)]
    results[experiment_type] = temp
    print '{0:30} {1:>5}  out of {2:>5}'.format(
      experiment_type,len(temp),len(experiment_type_dict[experiment_type]))
    if experiment_type == 'X-RAY DIFFRACTION':
      i = 0
      for x in temp:
        if data_dict2[x][1] < data_dict[x][2]:
          outstr = '{0:6} PROBE:{1:<5.1f} RM:{2:<7.1f} Sym:{3:<7.1f} Solvent:{4:<5.1f}'
          print outstr.format(x,data_dict[x][2],data_dict2[x][1],data_dict2[x][2],data_dict2[x][3])
          i += 1
          if i>15: break
      #for x in temp[50:55]:
        #print x,data_dict[x]
  print '='*85
  print 'The number of files with abs(clashscore - nb_clashscore) > {}: '.format(delta)
  print '-'*85
  for experiment_type in experiment_type_dict:
    temp = [x for x in experiment_type_dict[experiment_type] 
            if data_dict2.has_key(x) and data_dict.has_key(x) and 
            (abs(data_dict2[x][0]-data_dict[x][2]) > delta)]
    results[experiment_type] = temp
    print '{0:30} {1:>5}  out of {2:>5}'.format(
      experiment_type,len(temp),len(experiment_type_dict[experiment_type]))
  print '='*85
  print 'The number of files with abs(clashscore - nb_clashscore)/(clashscore + 0.001) > {}: '.format(ratio)
  print '-'*85
  for experiment_type in experiment_type_dict:
    temp = [x for x in experiment_type_dict[experiment_type] 
            if data_dict2.has_key(x) and data_dict.has_key(x) and
            (abs(data_dict2[x][0]-data_dict[x][2])/(data_dict[x][2] + 0.001)> ratio)]
    results[experiment_type] = temp
    print '{0:30} {1:>5}  out of {2:>5}'.format(
      experiment_type,len(temp),len(experiment_type_dict[experiment_type]))
  print '='*85
  print 'The number of files with abs(total_nb_clashscore - without_sym_nb_clashscore) > {}: '.format(delta)
  print '-'*85
  for experiment_type in experiment_type_dict:
    temp = [x for x in experiment_type_dict[experiment_type] 
            if data_dict2.has_key(x) and data_dict.has_key(x) and
            (abs(data_dict2[x][0]-data_dict2[x][1]) > delta)]
    results[experiment_type] = temp
    print '{0:30} {1:>5}  out of {2:>5}'.format(
      experiment_type,len(temp),len(experiment_type_dict[experiment_type]))
  print '='*85
  print 'The number of files with solvent-solvent clashscore > {}: '.format(delta_solvent)
  print '-'*85
  for experiment_type in experiment_type_dict:
    temp = [x for x in experiment_type_dict[experiment_type] 
            if data_dict2.has_key(x) and data_dict.has_key(x) and
            (data_dict2[x][3] > delta_solvent)]
    results[experiment_type] = temp
    print '{0:30} {1:>5}  out of {2:>5}'.format(
      experiment_type,len(temp),len(experiment_type_dict[experiment_type]))
    if experiment_type == 'X-RAY DIFFRACTION':
      i = 0
      for x in temp:
        outstr = '{0:6} PROBE:{1:<5.1f} RM:{2:<7.1f} Sym:{3:<7.1f} Solvent:{4:<5.1f}'
        print outstr.format(x,data_dict[x][2],data_dict2[x][1],data_dict2[x][2],data_dict2[x][3])
        i += 1
        if i>15: break
  print '='*85
  
  
  print 'Done...'
  
if __name__=='__main__':
  run()