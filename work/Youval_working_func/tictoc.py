'''
Matalab like tic() toc()
'''

import time

def tic():
    #Homemade version of matlab tic and toc functions
    global startTime_for_tictoc
    startTime_for_tictoc = time.time()

def toc(msg='',print_time=True):
    if 'startTime_for_tictoc' in globals():
        if print_time:
            outstr = '{0}: Elapsed time is: {1:.4f} seconds\n'.format(msg,time.time() - startTime_for_tictoc)
            print outstr
        else:
            outstr = '{0:.4f}'.format(time.time() - startTime_for_tictoc)
            return outstr
    else:
        print "Toc: start time not set"