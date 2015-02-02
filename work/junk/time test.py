from __future__ import division
import time

__author__ = 'Youval'


a = range(250*250)
b = {x:x for x in a}


t0 = time.time()
for i in xrange(10000000):
  a[10000]
  if a[200] > 0:
    a[100] += 1
print round(time.time()-t0,2)
t0 = time.time()
for i in xrange(10000000):
  b[10000]
  if b[200] > 0:
    b[100] += 1
print round(time.time()-t0,2)

t0 = time.time()
for i in xrange(10000000):
  a[10000]
  a[100] += 1 * (a[200] > 0)
print round(time.time()-t0,2)
t0 = time.time()
for i in xrange(10000000):
  b[10000]
  b[100] += 1 * (b[200] > 0)
print round(time.time()-t0,2)
