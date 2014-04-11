from __future__ import division
from scitbx.array_family import flex
import numpy as np



x = flex.random_double(12)
x = flex.vec3_double(x)
y = flex.random_double(6)*0
y = flex.vec3_double(y)

a = np.ones([10,3])
b = np.zeros([2,3])





# x.as_list()[6:10] = y.as_list()

print list(x[6:10])

x[6:10]=y


print list(x)


print 'done'