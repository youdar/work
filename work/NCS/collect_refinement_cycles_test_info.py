from __future__ import division
import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np

def run():
  test_data = open('translation_test_data.txt','r').read().splitlines()
  r_work_start = []
  r_work_end = []
  for rec in test_data:
    if rec.startswith('start r_factor:'):
      r_work_start.append(float(rec.split(':')[-1]))
    elif rec.startswith('  macro_cycle   9'):
      r_work_end.append(float(rec.split(' ')[10]))
  assert len(r_work_start)==len(r_work_end)

  s1 = 'R-work start mean: {0:.4f}   standard deviation: {1:0.4f}'
  s2 = 'Number of iteration that went down to R-work = 0 : {}'
  r_work_mean = np.mean(r_work_start)
  r_work_std = np.std(r_work_start)
  s1 = s1.format(r_work_mean, r_work_std)
  print s1
  # print s2.format(r_work_end.count(0))

  plot_histogram(data=r_work_end, n_bins=100, title_str=s1)

def plot_histogram(data,n_bins,title_str=''):
  """  Plot Histogram  """
  fig, ax = plt.subplots()

  # histogram our data with numpy
  n, bins, patches = ax.hist(data, bins=n_bins)

  # we need to normalize the data to 0..1 for the full
  # range of the colormap
  fracs = n.astype(float)/n.max()
  norm = colors.Normalize(fracs.min(), fracs.max())

  for thisfrac, thispatch in zip(fracs, patches):
      color = cm.jet(norm(thisfrac))
      thispatch.set_facecolor(color)

  plt.title(title_str)
  plt.xlabel('R-work')
  plt.ylim((0,1.1*n.max()))
  # plt.xlim((0,max(0.1,bins.max())))
  plt.show()



if __name__=='__main__':
  run()