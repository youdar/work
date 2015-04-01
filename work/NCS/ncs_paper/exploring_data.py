from __future__ import division
from pandas.tools.plotting import scatter_matrix
import matplotlib.pyplot as plt
from libtbx.utils import Sorry
import scipy.optimize as op
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
    if not os.path.isfile(fn):
      raise Sorry('Make sure ncs_paper_data.csv is in %s'%c.ncs_dir)
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


  def plot_final_values(self,var_name):
    """
    Args:
      var_name (str): the name of the variable we want to make a plot for
    """
    print '-'*50
    print 'Making plot for:',var_name
    plt.close("all")
    rfree_list = ['{} : {}'.format(var_name,x) for x in self.refine_test_names]
    rfreefinal = self.df[rfree_list]
    rfinal = rfreefinal.dropna(axis=0)
    print 'Total number of files (including with no results:',len(rfreefinal)
    print 'Total number of files:',len(rfinal)
    # rfinal = rfinal.sort(['r-free final : no ncs'], ascending=True)
    rfinal = rfinal[rfinal > 0].dropna()
    print 'Number of plotted files: ',len(rfinal)
    print rfinal.head(5)
    f = plt.figure(figsize=(10,6))
    # Plot r-free differene between methods
    x = range(len(rfinal))
    # sort values in by the 'no ncs' order
    values = (rfinal.values).transpose()
    values = values[:,np.argsort(values[0,:])]
    plot_types = ['.-b','.y','.k']
    opacity = [1,0.5,0.5]
    assert rfree_list[0].split(' : ')[1] == 'no ncs'
    for refine_t in rfree_list:
        i = rfreefinal.columns.get_loc(refine_t)
        y = values[i]
        plt.plot(x,y,plot_types[i],alpha=opacity[i])
    # add text to plot
    txt1 = 'Number of structures: {}'.format(len(rfinal))
    plt.annotate(txt1,xy=(0.7*max(x),min(y)))
    if var_name == 'final clashscore' : y_lable = 'Clashscore'
    if var_name == 'r-free final' : y_lable = 'R-Free Final'
    plt.xlabel('PDB structures')
    plt.ylabel(y_lable)
    plt.legend(self.refine_test_names,fontsize=14,loc=2)
    fig_name = y_lable + '.png'
    c = collect_ncs_files.ncs_paper_data_collection()
    plt.savefig(os.path.join(c.figures_dir,fig_name),ext="png",dpi=600)
    plt.show()


  def save_pdb_ids_with_ncs_issues(self,df):
    """
    When running refinement with the default value for
    refinement.ncs.excessive_distance_limit

    Many files had issues with Excessive distances to NCS averages,
    and refinement did not complete.

    We changed the refinement option to:
    Excessive distances to NCS averages=None

    But saved the PDB IDs with issues in:
    files_with_excessive_distance_limit_issues.txt
    """
    g = df['r-free final : no ncs'] > 0
    df2 = df[[
      'pdb id',
      'r-free final : cartesian ncs restraints',
      'r-free final : torsion ncs restraints']][g]
    g2 = (df2['r-free final : cartesian ncs restraints'] == 0)
    g2 |= (df2['r-free final : torsion ncs restraints'] == 0)
    pdb_ids = df2['pdb id'][g2].values
    # Save file list
    c = collect_ncs_files.ncs_paper_data_collection()
    pdb_ids = '\n'.join(pdb_ids)
    fn = os.path.join(c.ncs_dir,'files_with_excessive_distance_limit_issues.txt')
    open(fn,'w').write(pdb_ids)

  def find_outliers(self,df,var_name,test_name):
    """
    find cases where the value of 'var_name : test_name' is larger than
    'var_name : no ncs'

    Args:
      df : data frame
      var_name (str): the name of the variable we want to make a plot for
      test_name (str): refinement test name
    """
    # create the names of the columns to use
    tests = ['no ncs',test_name]
    col_list = ['{} : {}'.format(var_name,x) for x in tests]
    col_list.insert(0,'pdb id')
    print col_list
    # Get those columns
    f = self.get_clean_data_frame(df)
    f = f[col_list]
    f = f.dropna()
    g = f[col_list[1]]<f[col_list[2]]
    f = f[g]
    f['delta'] = f[col_list[2]]  - f[col_list[1]]
    f = f.sort(['delta'],ascending=False)
    print 'number of files with worse results:',sum(g)
    print f.head(6)
    print f.describe()
    print 'done'


  def get_clean_data_frame(self,df):
    """ remove rows with NA and zeros"""
    g = df['r-free final : no ncs'] > 0
    g &= (df['r-free final : cartesian ncs restraints'] > 0)
    g &= (df['r-free final : torsion ncs restraints'] > 0)
    return df[g]

  def find_when_to_use_ncs(self,df,var_name,test_name):
    """
    Learn from data what are the conditions where NCS refinement will be most
    useful

    Args:
      df : data frame
      var_name (str): the name of the variable we want to make a plot for
      test_name (str): refinement test name

    Returns:

    """
    # get the results:
    #   y = 1 : NCS improved results (smaller value)
    #   y = 0 : NCS did not improve
    tests = ['no ncs',test_name]
    col_list = ['{} : {}'.format(var_name,x) for x in tests]
    col_list.insert(0,'pdb id')

    # drop non numerical columns
    # todo: check why the following have missing values (look at mtz files)
    '''
    3asn

    '''
    drop_list = ['pdb id','year','r-free header',
                 'r-work header','master only',
                 'experiment',
                 'r-free init : torsion ncs restraints',
                 'r-free init : cartesian ncs restraints']

    remove_n = [
      'refinement time','r-work final',
      'cbeta deviations','cbeta final','rotamer outliers','rama outliers',
      'rotamer final','r-work init']
    # replace column names with shorter names
    new_names_1=[
      'Copies','Groups','Res.','Comp.','Solv.',
      'ASU Atom','p/d ncs','p/d asu',
      'R-Free init', 'R-Free final',
      'R final Cartesian',
      'R final Torsion']
    old_names_1=[
      'n copies', 'n groups', 'resolution', 'completeness', 'solvent fraction',
      'atoms in asu', 'p-to-d ratio ncs', 'p-to-d ratio asu',
      'r-free init : no ncs', 'r-free final : no ncs',
      'r-free final : cartesian ncs restraints',
      'r-free final : torsion ncs restraints']
    remove_n_1 = remove_n + ['final clashscore','all-atom clashscore']

    n = 100
    sampled_f = self.get_partial_df(
      df,
      drop_list=drop_list,
      new_names=new_names_1,
      old_names=old_names_1,
      remove_n=remove_n_1,
      n=n)
    scatterplot_matrix(sampled_f,fig_name='Grid_R_free.png')

    test_df = pd.DataFrame(np.random.randn(1000,4),columns=['a','b','c','d'])
    # scatterplot_matrix(test_df)

    # plot 2
    new_names_2=[
      'Copies','Groups','Res.','Comp.','Solv.',
      'ASU Atom','p/d ncs','p/d asu',
      'R-Free final',
      'clashsocre',
      'C.S. Cartesian',
      'C.S. Torsion']
    old_names_2=[
      'n copies', 'n groups', 'resolution', 'completeness', 'solvent fraction',
      'atoms in asu', 'p-to-d ratio ncs', 'p-to-d ratio asu',
      'r-free final : no ncs',
      'final clashscore : no ncs',
      'final clashscore : cartesian ncs restraints',
      'final clashscore : torsion ncs restraints']
    drop_list_2 = drop_list + ['r-free final : cartesian ncs restraints']
    drop_list_2.append('r-free final : torsion ncs restraints')
    drop_list_2.append('r-free init : no ncs')
    remove_n_2 = remove_n + ['all-atom clashscore']

    sampled_f = self.get_partial_df(
      df,
      drop_list=drop_list_2,
      new_names=new_names_2,
      old_names=old_names_2,
      remove_n=remove_n_2,
      n=n)
    scatterplot_matrix(sampled_f,fig_name='Grid_clashscore.png')

    # g = f[col_list[1]]>f[col_list[2]]
    # # remove the answer column from data
    # f.drop(col_list[2],axis=1,inplace=True)
    # y = np.array(g*1)
    # y = y.reshape((len(g),1))
    print 'Done'

  def get_partial_df(self,df,drop_list,new_names,old_names,remove_n,n):
    """ remove item from dataframe """
    drop_list = list(drop_list)
    refine_n = collect_ncs_files.get_refine_test_names()
    f = self.get_clean_data_frame(df)
    f = f.dropna(axis=1)
    for tst in refine_n:
      remove_list = ['{} : {}'.format(x,tst) for x in remove_n]
      drop_list.extend(remove_list)
    for col in drop_list:
      if col in f.columns:
        f.drop(col,axis=1,inplace=True)
    rows = np.random.choice(f.index.values, n)
    f = f.ix[rows]
    assert len(new_names) == len(old_names)
    assert list(f.columns) == old_names
    f.columns = new_names
    return f

def scatterplot_matrix(df,fig_name=''):
  g = sns.PairGrid(
    df,
    hue="p/d ncs",
    size=1.2,
    # aspect=1.6,
    dropna=True)
  # fig.square_grid = True
  g.map_diag(plt.hist)
  g.map_offdiag(plt.scatter)
  # reduce the number of values in the legend
  legend_keys = sorted(g._legend_data)
  i = len(legend_keys) - 1
  d = i//10 + 1
  keys = []
  while i>=0:
    keys.append(legend_keys[i])
    i -= d
  legend = {k:g._legend_data[k] for k in keys}
  g.add_legend(title='Data/Param NCS',legend_data=legend)
  g.set(ylim=(0, None))
  g.set(xlim=(0, None))
  # set the space between subplot
  g.fig.subplots_adjust(wspace=0.02,hspace=0.02)
  # set the number of ticks in each subplot
  allticks = g.fig.get_axes()
  for ticks in allticks:
    # reduce the number of ticks
    tx = ticks.get_xticks()
    ty = ticks.get_yticks()
    if len(tx) > 3:
      mx = len(tx)//2
      ticks.set_xticks([tx[0],tx[mx],tx[-2]])
    if len(ty) > 3:
      my = len(ty)//2
      ticks.set_yticks([ty[0],ty[my],ty[-2]])
  # fig.map(plt.scatter)
  if fig_name:
    c = collect_ncs_files.ncs_paper_data_collection()
    plt.savefig(os.path.join(c.figures_dir,fig_name),ext="png",dpi=300)
  # plt.show()

def learn_relations_param_to_good_ncs_effect(data):
  """
  Use logistic regression to evaluate when will the NCS refinement results
  will not be as good as those without NCS

  Args:
    data : training data m-by-n matrix, where the last
      column is the 1: < no-ncs, 0: > no-ncs
  """
  # set y as the last column, the "answers"
  y = data[:,-1]
  # set x as the data
  X = data[:,:-1]
  # Normalize X parameters

  m,n = X.shape
  y = y.reshape((m,1))
  # add a column with the value 1 at the start of the parameters
  X = np.append(np.ones([m,1]),X,1)
  # Initialize fitting parameters
  initial_theta = np.zeros([n + 1, 1])
  # Compute and display initial cost and gradient
  cost = cost_function(initial_theta, X, y)
  grad = gradient(initial_theta, X, y)
  print 'Initial cost:',cost
  # Call minimizer
  options = {'full_output': True, 'maxiter': 400}
  myargs = (X, y)
  out = op.fmin_bfgs(
    cost_function,
    initial_theta,
    args=myargs,
    fprime=gradient,
    **options)
  # collect results
  optimal_theta,cost,grad_at_min,inv_hessian,fun_calls,grad_calls,flags = out
  print 'Cost after minimization:',cost
  print 'Parameters:',optimal_theta
  # return predictor with optimal theta
  return predict_func(optimal_theta)

def cost_function(theta, X, y):
  """
  Compute cost and gradient for logistic regression

  Args:
    theta (array) : n-by-1 array of current parameters values
    X (2D array): training data
    y (array): training answers

  Returns:
    cost (float): costFunction(theta, X, y) computes the cost of using theta
  """
  #Initialize some useful values
  m = len(y) # number of training examples
  h = sigmoid(X.dot(theta))
  cost = -(y.transpose().dot(np.log(h))+(1-y).transpose().dot(np.log(1-h)))/m
  grad = X.transpose().dot((h - y))/m
  return cost[0,0], grad

def gradient(theta,X,y):
  """
  Return the gradient of a logistic regression
  grad (array) : the gradient of the cost with the respective theta parameters
  """
  h = sigmoid(X.dot(theta))
  m,n = X.shape
  h = h.reshape((m,1))
  y = y.reshape((m,1))
  grad = X.transpose().dot((h - y))/m
  return grad

def sigmoid(x):
  """ returns a Sigmoid of an array x """
  return 1/(1 + np.exp(-x))

def predict_func(theta):
  def predict(x):
    h = sigmoid(x.dot(theta))
    return np.round(h,2)
  return predict


def run():
  explore = Explore_data()
  df = explore.get_data_frame()
  #
  explore.find_when_to_use_ncs(df,'r-free final','cartesian ncs restraints')
  # explore.plot_delta_r_values(['r-work final','r-free final'])
  # explore.find_outliers(df,'r-free final','cartesian ncs restraints')
  print '+'*50
  explore.find_outliers(df,'final clashscore','cartesian ncs restraints')
  print '+'*50
  # #
  # explore.plot_final_values('r-free final')
  # explore.plot_final_values('final clashscore')
  # explore.save_pdb_ids_with_ncs_issues(df)
  print 'Done'


if __name__ == '__main__':
  run()









