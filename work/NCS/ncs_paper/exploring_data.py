from __future__ import division
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
    tests = ['no ncs',test_name]
    col_list = ['{} : {}'.format(var_name,x) for x in tests]
    col_list.insert(0,'pdb id')
    print col_list

    f = self.get_clean_data_frame(df)
    f = f[col_list]
    f = f.dropna()
    g = f[col_list[1]]<f[col_list[2]]
    f = f[g]
    f['delta'] = f[col_list[2]]  - f[col_list[1]]
    f = f.sort(['delta'],ascending=False)
    print len(f)
    print f.head(6)

    print 'done'


  def get_clean_data_frame(self,df):
    """ remove rows with NA and zeros"""
    g = df['r-free final : no ncs'] > 0
    g &= (df['r-free final : cartesian ncs restraints'] > 0)
    g &= (df['r-free final : torsion ncs restraints'] > 0)
    return df[g]


def learn_relations_param_to_bad_refine(data):
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
  # explore.plot_delta_r_values(['r-work final','r-free final'])
  explore.find_outliers(df,'r-free final','cartesian ncs restraints')
  print '+'*50
  explore.find_outliers(df,'final clashscore','cartesian ncs restraints')
  print '+'*50
  #
  explore.plot_final_values('r-free final')
  explore.plot_final_values('final clashscore')
  # explore.save_pdb_ids_with_ncs_issues(df)
  print 'Done'


if __name__ == '__main__':
  run()









