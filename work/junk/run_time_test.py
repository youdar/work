from __future__ import division
from scitbx import matrix
import time

def g(x,y):
  return 0.5*(1 + x + y)*(x + y) + y

def f(r):
  (r11,r12,r13,r21,r22,r23,r31,r32,r33) = (10*r.transpose()).round(0).elems
  return str(int(g(r11,g(r22,r33))))

def run():
  for i in range(120):
    for j in range(i,120):
      r = matrix.sqr([0.1235,1.0,1.0,0.2,0.5,0.6,0.7,0.8,0.9])
      print f(r)

if __name__ == "__main__":
  t0=time.time()
  run()
  print "Time: %6.4f"%(time.time()-t0)
  print "OK"

