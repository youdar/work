import numpy as nm
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from  iotbx.pdb.multimer_reconstruction import *
import os


def run():
    '''
    
    
    '''
    os.chdir('/net/cci/youval/Work/')
    print(os.getcwd())
    v = multimer('1S58.pdb','ba')
    
    

    v.write()
    print os.getcwd() + '/' + v.pdb_output_file_name
    print 'ok'
    
    

if __name__=="__main__":
    run()
    
    
