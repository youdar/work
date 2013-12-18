from __future__ import division
import cPickle as pickle
import os

'''
Collecting pdb files that have very different clashscore and nb_clashscore
So we can look at the cause and evaluate the difference in scoring methods
'''


  
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

  return clashscore_data_dict, experiment_type_dict

def run():
  clashscore_data_dict, experiment_type_dict = get_data()
  # Collect files where the abs(clashscore - nb_clashscore) > 10
  results = {}
  delta = 30
  delta_solvent = 10
  ratio = 0.3
  year_limit = 2000
  #
  print 'nb_simple is nonbonded clashscore without sym op and solvent-solvent clashes'
  print 'The number of files with abs(clashscore (O_vdw 1.52) - nb_simple ) > {}: '.format(delta)
  print '-'*85
  for experiment_type in experiment_type_dict:
    temp = [x for x in experiment_type_dict[experiment_type] 
            if clashscore_data_dict.has_key(x) and 
            clashscore_data_dict[x][6]>year_limit and
            (abs(clashscore_data_dict[x][1]-clashscore_data_dict[x][5]) > delta)]
    results[experiment_type] = temp
    print '{0:30} {1:>5}  out of {2:>5}'.format(
      experiment_type,len(temp),len(experiment_type_dict[experiment_type]))
    #
    print_info(experiment_type='X-RAY DIFFRACTION',
               temp=temp,
               clashscore_data_dict=clashscore_data_dict,
               n=15)
  print '='*85
  #
  print 'The number of files with abs(clashscore (O_vdw 1.52)- Total nb_clashscore) > {}: '.format(delta)
  print '-'*85
  for experiment_type in experiment_type_dict:
    temp = [x for x in experiment_type_dict[experiment_type] 
            if clashscore_data_dict.has_key(x) and clashscore_data_dict[x][6]>year_limit and 
            (abs(clashscore_data_dict[x][0]-clashscore_data_dict[x][5]) > delta)]
    results[experiment_type] = temp
    print '{0:30} {1:>5}  out of {2:>5}'.format(
      experiment_type,len(temp),len(experiment_type_dict[experiment_type]))
  print '='*85
  #
  print 'The number of files with abs(clashscore (O_vdw 1.52) - nb_simple)/(clashscore + 0.001) > {}: '.format(ratio)
  print '-'*85
  for experiment_type in experiment_type_dict:
    temp = [x for x in experiment_type_dict[experiment_type] 
            if clashscore_data_dict.has_key(x) and clashscore_data_dict[x][6]>year_limit and
            (abs(clashscore_data_dict[x][1]-clashscore_data_dict[x][5])/(clashscore_data_dict[x][5] + 0.001)> ratio)]
    results[experiment_type] = temp
    print '{0:30} {1:>5}  out of {2:>5}'.format(
      experiment_type,len(temp),len(experiment_type_dict[experiment_type]))
  print '='*85
  #
  print 'The number of files with abs(total_nb_clashscore - nb_simple) > {}: '.format(delta)
  print '-'*85
  for experiment_type in experiment_type_dict:
    temp = [x for x in experiment_type_dict[experiment_type] 
            if clashscore_data_dict.has_key(x) and clashscore_data_dict[x][6]>year_limit and
            (abs(clashscore_data_dict[x][0]-clashscore_data_dict[x][1]) > delta)]
    results[experiment_type] = temp
    print '{0:30} {1:>5}  out of {2:>5}'.format(
      experiment_type,len(temp),len(experiment_type_dict[experiment_type]))
  print '='*85
  #
  print 'The number of files with solvent-solvent clashscore > {}: '.format(delta_solvent)
  print '-'*85
  for experiment_type in experiment_type_dict:
    temp = [x for x in experiment_type_dict[experiment_type] 
            if clashscore_data_dict.has_key(x) and clashscore_data_dict[x][6]>year_limit and
            (clashscore_data_dict[x][3] > delta_solvent)]
    results[experiment_type] = temp
    print '{0:30} {1:>5}  out of {2:>5}'.format(
      experiment_type,len(temp),len(experiment_type_dict[experiment_type]))
    #
    print_info(experiment_type='X-RAY DIFFRACTION',
               temp=temp,
               clashscore_data_dict=clashscore_data_dict,
               n=15)
  print '='*85
  
  print 'Done...'
  
def print_info(experiment_type,temp,clashscore_data_dict,n=15):
  if experiment_type == 'X-RAY DIFFRACTION':
      i = 0
      for x in temp:
        if clashscore_data_dict[x][1] < clashscore_data_dict[x][2]:
          outstr = '{0:6} PROBE (O_vdw 1.52):{1:<5.1f} RM (simple):{2:<7.1f} Sym:{3:<7.1f} Solvent:{4:<5.1f}'
          print outstr.format(x,
                              clashscore_data_dict[x][5],
                              clashscore_data_dict[x][1],
                              clashscore_data_dict[x][2],
                              clashscore_data_dict[x][3])
          i += 1
          if i>n: break

  
if __name__=='__main__':
  run()