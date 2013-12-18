import os
import cPickle as pickle

'''
combine the data from
clashscore_compare_reduce_12_6_2013_dict     dict[pdb_file_name] = [total_nb_clashscore,without_sym_nb_clashscore,clashscore_probe]
clashscore_compare_reduce_12_11_2013_dict    dict[pdb_file_name] = [total_nb_clashscore,without_sym_nb_clashscore,clashscore_probe]
clashscore_solvent_reduce_12_13_2013_dict    dict[pdb_file_name] = [clashscore_all_clashes,clashscore_simple,clashscore_only_sym_op,clashscore_solvent_solvent]
file_to_year_dict

using
experiment_type_to_files_dict

to one dictionary, clashscore_data_dict
clashscore_data_dict[pdb_file] = [clashscore_all_clashes,
                                  clashscore_simple,
                                  clashscore_only_sym_op,
                                  clashscore_solvent_solvent
                                  probe_clashscore_Ovdw_14,
                                  probe_clashscore_Ovdw_152,
                                  pdb_year]
'''

def run():
  # set path
  datapath = os.path.realpath('c:\Phenix\Dev\Work\work\Clashes\Data')
  # Read data
  data_dict_file_Ovdw_140 = 'clashscore_compare_reduce_12_6_2013_dict' 		# Probe O vdw is 1.4
  data_dict_Ovdw_140 = pickle.load(open(os.path.join(datapath,data_dict_file_Ovdw_140),'r'))
  
  data_dict_file_Ovdw_152 = 'clashscore_compare_reduce_12_11_2013_dict' 	# Probe O vdw is 1.52
  data_dict_Ovdw_152 = pickle.load(open(os.path.join(datapath,data_dict_file_Ovdw_152),'r'))
  
  data_dict_file_solvent = 'clashscore_solvent_reduce_12_13_2013_dict' 		# with solvent-solvent results
  data_dict_solvent = pickle.load(open(os.path.join(datapath,data_dict_file_solvent),'r'))
  
  experiment_dict_file = 'experiment_type_to_files_dict'		# source for experiment_type_dict
  experiment_type_dict = pickle.load(open(os.path.join(datapath,experiment_dict_file),'r'))
  
  file_to_year_dict = pickle.load(open(os.path.join(datapath,'file_to_year_dict'),'r'))
  #
  # files that are not yet in the index
  missing_files = []
  
  clashscore_data_dict = {}
  for experiment_type in experiment_type_dict:
    for pdb_file in experiment_type_dict[experiment_type]:
      if file_to_year_dict.has_key(pdb_file):
        if data_dict_Ovdw_140.has_key(pdb_file) and \
           data_dict_Ovdw_152.has_key(pdb_file) and \
           data_dict_solvent.has_key(pdb_file):
          s = []
          s = data_dict_solvent[pdb_file]
          s.append(data_dict_Ovdw_140[pdb_file][-1])
          s.append(data_dict_Ovdw_152[pdb_file][-1])
          s .append(file_to_year_dict[pdb_file])
          clashscore_data_dict[pdb_file] = s
        #else:
          #missing_files.append(pdb_file)
          
  print len(missing_files)     
  #pickle.dump(missing_files,open(os.path.join(datapath,'missing_files'),'w'))
  pickle.dump(clashscore_data_dict, open(os.path.join(datapath,'clashscore_data_dict'),'w'))
  print 'Done...'
  
if __name__=='__main__':
  run()
  