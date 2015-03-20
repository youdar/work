from __future__ import division
import matplotlib.pyplot as plt
import collect_ncs_files
import seaborn as sns
import pandas as pd
import numpy as np
import sys
import os


class Explore_data(object):

  def __init__(self):
    self.df = None
    headers, table_pos_map = collect_ncs_files.table_headers()
    self.headers = headers
    self.table_pos_map = table_pos_map
    self.refine_test_names = collect_ncs_files.get_refine_test_names()
    # Del: just for testing
    self.refine_test_names = [
      'no ncs','cartesian ncs restraints','torsion ncs restraints']
    # experiment and plot types
    self.plot_types = ['r-free final','final clashscore']
    self.plot_deltas = ['r-work final','r-free final']

  def get_data_frame(self):
    """ Create pandas data frame """
    c = collect_ncs_files.ncs_paper_data_collection()
    fn = os.path.join(c.ncs_dir,'ncs_paper_data.csv')
    print os.path.isfile(fn)
    df = pd.pandas.DataFrame.from_csv(fn,index_col=False)
    self.df = df
    return df

  def plot_clashscore(self):
    pass

  def plot_delta_r_values(self,vals):
    """
    Args:
      vals (list): list of length 2, with the columns to be compared
    """
    print 'Plotting difference scatter plot for: ({} , {})'.format(*vals)
    df = self.df
    rfree_list = ['{} : {}'.format(vals[1],x) for x in self.refine_test_names]
    rwork_list = ['{} : {}'.format(vals[0],x) for x in self.refine_test_names]
    # look at r-free final
    rfreefinal = df[rfree_list]
    # look at r-work final
    rworkfinal = df[rwork_list]
    # look at the difference
    x = rworkfinal.values
    y = rfreefinal.values
    # get values for r_work - r_free
    rdiff = pd.concat(
      [df['pdb id'],pd.DataFrame(x-y,columns=self.refine_test_names)],axis=1)
    # relative values
    c = 'r-free final : no ncs'
    i = rfreefinal.columns.get_loc('r-free final : no ncs')
    # remove all PDBs without refinement results
    rdiff = rdiff.dropna(axis=0)
    print rdiff.head(5)
    print 'Number of files with all results:',len(rdiff)
    # Scatter Plot r-work , r-free difference
    values = (rdiff.values).transpose()
    y = values[1]
    plot_types = ['.b','.y','.k']
    opacity = [1,0.5,0.5]
    for i in range(1,len(self.refine_test_names)):
        print self.refine_test_names[i],i,plot_types[i],opacity[i]
        x = values[i]
        plt.plot(x,y,plot_types[i],alpha=opacity[i])

    min_x = min(y)
    max_x = max(y)
    print min_x,max_x
    plt.plot([min_x,max_x],[min_x,max_x],'-b')
    plt.xlim([min_x,max_x])
    plt.ylim([min_x,max_x])
    plt.show()

    # plot
    # sns.tsplot(df, err_style="unit_points", color="mediumpurple")
    # df.plot(kind='',figsize=(15,5))
    # add text to plot
    # plt.annotate()


  def plot_r_free_finals(self,var_name):
    """
    Args:
      var_name (str): the name of the variable we want to make a plot for
    """
    print 'Making plot for:',var_name
    rfree_list = ['{} : {}'.format(var_name,x) for x in self.refine_test_names]
    rfreefinal = self.df[rfree_list]
    rfinal = rfreefinal.dropna(axis=0)
    # rfinal = rfinal.sort(['r-free final : no ncs'], ascending=True)
    print rfinal.head(5)

    # Plot r-free differene between methods
    x = range(len(rfinal))
    values = (rfinal.values).transpose()
    values = values[:,np.argsort(values[0,:])]
    plot_types = ['.-b','.y','.k']
    opacity = [1,0.5,0.5]
    for refine_t in rfree_list:
        i = rfreefinal.columns.get_loc(refine_t)
        y = values[i]
        plt.plot(x,y,plot_types[i],alpha=opacity[i])
    # add text to plot
    # plt.annotate()
    plt.show()


def run():
  explore = Explore_data()
  df = explore.get_data_frame()
  # print df.head(3)
  a = df[['pdb id','n groups']]
  explore.plot_delta_r_values(['r-work final','r-free final'])
  # explore.plot_r_free_finals('r-free final')


  print 'Done'


if __name__ == '__main__':
  run()









