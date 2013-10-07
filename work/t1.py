#from multiprocessing import Process, Lock

#def f(l, i):
    #l.acquire()
    #print 'hello world', i
    #l.release()

#if __name__ == '__main__':
    #lock = Lock()

    #for num in range(10):
        #Process(target=f, args=(lock, num)).start()

from __future__ import division
from scitbx import matrix
from  iotbx.pdb.multimer_reconstruction import multimer
import os

def run(file_name):
    cau_multimer_data = multimer(file_name,'cau')
    print 'done'

if __name__=='__main__':
    os.chdir('/net/cci-filer2/raid1/home/youval/Work')
    file_name = '01.pdb'
    if os.path.isfile(file_name):
        run(file_name)
    else:
        print 'No Such file'