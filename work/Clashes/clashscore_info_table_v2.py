from __future__ import division
import cPickle as pickle
import os

'''
Collecting pdb files that have very different clashscore and nb_clashscore
So we can look at the cause and evaluate the difference in scoring methods
'''

########################################################################
class data_table_plot(object):
  """"""
  #----------------------------------------------------------------------
  def __init__(self):
    """Constructor"""
    self.delta = 10
    self.delta_solvent = 10
    self.ratio = 0.3
    self.year_limit = 2000
    self.smaple_size = 15
    print 'nb_simple is nonbonded clashscore without sym op and solvent-solvent clashes'
    print 'data for pdb file with no unknown pairs data, from the year {} and on'.format(self.year_limit)

    

  def get_data(self):
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
    
    self.experiment_type_dict = experiment_type_dict
    self.clashscore_data_dict = clashscore_data_dict

  
  
  def print_table(self,title,i_rec,j_rec,test_type,test_val,print_sample_experiment_type=[]):
    '''
    rec number:
    0: clashscore_all_clashes
    1: clashscore_simple,
    2: clashscore_only_sym_op,
    3: clashscore_solvent_solvent
    4: probe_clashscore_Ovdw_14,
    5: probe_clashscore_Ovdw_152,
    6: pdb_year
    '''
    results = {}
    outstr = '{0:6} PROBE (1.52):{1:<5.1f} (1.4):{2:<5.1f} RM (simple):{3:<7.1f} Sym:{4:<7.1f} Solvent:{5:<5.1f}  : {6}'
    print '='*85
    print title
    print '-'*85
    for experiment_type in self.experiment_type_dict:
      records = [x for x in self.experiment_type_dict[experiment_type] 
              if self.clashscore_data_dict.has_key(x) and self.clashscore_data_dict[x][6]>=self.year_limit]
      temp = []
      i = 0
      for x in self.experiment_type_dict[experiment_type]:
        if self.clashscore_data_dict.has_key(x) and self.clashscore_data_dict[x][6]>=self.year_limit and \
           self.test(val1=self.clashscore_data_dict[x][i_rec],
                     val2=self.clashscore_data_dict[x][j_rec],
                     test_type=test_type,test_val=test_val):
          temp.append(x)
          if (i < self.smaple_size) and experiment_type in print_sample_experiment_type:
            i += 1
            print outstr.format(x,
                              self.clashscore_data_dict[x][5],
                              self.clashscore_data_dict[x][4],
                              self.clashscore_data_dict[x][1],
                              self.clashscore_data_dict[x][2],
                              self.clashscore_data_dict[x][3],
                              experiment_type) 
          
        
      results[experiment_type] = temp
      print '{0:30} {1:>5}  out of {2:>5}'.format(experiment_type,len(temp),len(records))

        
  def test(self,val1,val2,test_type,test_val):
    if test_type == 'delta':
      return abs(val1-val2) > test_val
    elif test_type == 'signed_delta':
      return val1-val2 > test_val
    elif test_type == 'ratio':
      return abs(val1-val2)/(val2 + 0.001) > test_val
    elif test_type == 'None':
      return val1 > test_val
    else:
      raise TypeError,'Wrong test type'

def run():
  '''
    rec number:
    0: clashscore_all_clashes
    1: clashscore_simple,
    2: clashscore_only_sym_op,
    3: clashscore_solvent_solvent
    4: probe_clashscore_Ovdw_14,
    5: probe_clashscore_Ovdw_152,
    6: pdb_year
    '''
  table = data_table_plot()
  table.get_data()
  #
  title = 'The number of files with (clashscore (O_vdw 1.4) - total nb_clashscore) > {}: '.format(table.delta)
  table.print_table(title=title, test_type='signed_delta',i_rec=4,j_rec=0,
                    test_val=table.delta,print_sample_experiment_type=['X-RAY DIFFRACTION'])
  #
  title = 'The number of files with (clashscore (O_vdw 1.4) - nb_simple) > {}: '.format(table.delta)  
  table.print_table(title=title, test_type='signed_delta',i_rec=4,j_rec=1,
                    test_val=table.delta,print_sample_experiment_type=['X-RAY DIFFRACTION'])
  #
  title = 'The number of files with (nb_simple - clashscore (O_vdw 1.52)) > {}: '.format(table.delta)
  table.print_table(title=title, test_type='signed_delta',i_rec=1,j_rec=5,
                    test_val=table.delta,print_sample_experiment_type=['X-RAY DIFFRACTION'])
  #
  title =  'The number of files with abs(clashscore (O_vdw 1.52)- Total nb_clashscore) > {}: '.format(table.delta)
  table.print_table(title=title, test_type='delta',i_rec=0,j_rec=5,test_val=table.delta)
  #
  title = 'The number of files with abs(clashscore (O_vdw 1.52) - nb_simple)/(clashscore + 0.001) > {}: '.format(table.ratio)
  table.print_table(title=title, test_type='ratio',i_rec=1,j_rec=5,test_val=table.ratio)
  #
  title = 'The number of files with abs(total_nb_clashscore - nb_simple) > {}: '.format(table.delta)
  table.print_table(title=title, test_type='delta',i_rec=0,j_rec=1,test_val=table.delta)
  #
  title = 'The number of files with abs(clashscore (O_vdw 1.52) - clashscore (O_vdw 1.4)) > {}: '.format(table.delta)
  table.print_table(title=title, test_type='delta',i_rec=5,j_rec=4,test_val=table.delta)
  #
  title = 'The number of files with solvent-solvent clashscore > {}: '.format(table.delta_solvent)
  table.print_table(title=title, test_type='None',i_rec=3,j_rec=3,
                    test_val=table.delta_solvent,print_sample_experiment_type=[])
  #
  title = 'The number of files with symetry clashscore > {}: '.format(table.delta_solvent)
  table.print_table(title=title, test_type='None',i_rec=2,j_rec=2,
                    test_val=table.delta_solvent,print_sample_experiment_type=[])
  print '='*85
  print 'Done...'
  


  
if __name__=='__main__':
  run()