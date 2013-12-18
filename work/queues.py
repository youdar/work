from multiprocessing import Process, Queue, Lock
import time
from Just_wait import *
from tictoc import *



def reader(queue):
    while True:
        msg = queue.get()
        if (msg == 'DONE'):
            break

def reader2(queue, lock):
    while True:
        x = queue.get()
        time.sleep(.1)
        lock.acquire()
        print 'out:',x
        lock.release()
        if x == 'DONE':
            break


def writer(count, queue):
    for ii in xrange(count):
        queue.put(ii)
    queue.put('DONE')

if __name__=='__main__':
    lock = Lock()

    queue = Queue()
    reader_p = Process(target=reader2, args=(queue,lock))
    reader_p.daemon = True
    reader_p.start()	# launch the reader process
    tic()
    _start = time.time()	# the _ is just to symbolize previuos value
    for count in xrange(10**1):
        print count
        queue.put(count)
    queue.put('DONE')
    toc()

    reader_p.join()		# wait for the reader to finish
    toc()

    #for count in [10**4, 10**5, 10**6]:
        #queue = Queue()
        #reader_p = Process(target=reader, args=((queue),))
        ##
        #reader_p.daemon = True
        #reader_p.start()	# launch the reader process

        #_start = time.time()	# the _ is just to symbolize previuos value
        #writer(count, queue)	# send a lot of stuff to reader()
        #reader_p.join()		# wait for the reader to finish
        #print 'Sending %s numbers to Queue() took %s seconds' % (count,
                                                                 #(time.time() - _start))