#!/usr/bin/env ipython
#  Example 5 230 Notes  Lagrange points 1, 2, and 3

import numpy as np
import pylab
from scipy.optimize import brentq

mm = 7.3477e22  #  mass moon
me = 5.9736e24  #  mass Earth
r12 = 384400.0  #  distance between the Earth and the moon
pi2 = mm / (me + mm)


def f(xi):
    return ((1.0 - pi2) / np.absolute(xi + pi2) ** 3.0 * (xi + pi2) +
            pi2 / np.absolute(xi + pi2 - 1.0) ** 3.0 * (xi + pi2 - 1.0) - xi)


xi1 = brentq(f, -1.1, -1.0)
xi2 = brentq(f, 0.8, 0.9)
xi3 = brentq(f, 1.1, 1.2)

print "The zeros of the polynomial are", xi1, ",", xi2, ",", "and", xi3

print "L1 =", xi2 * r12
print "L2 =", xi3 * r12
print "L3 =", xi1 * r12

xi = np.linspace(-2.0, 2.0, 50000.0)

#  removing asymptotes
fxi = f(xi)
ext = [np.argmin(fxi), np.argmax(fxi)]
fxi[ext] = np.nan

fig = pylab.figure()
ax = fig.add_subplot(111)
ax.plot(xi, fxi, 'r')
pylab.grid()
pylab.xlim((-1.5, 1.5))
pylab.ylim((-1.0, 1.0))
pylab.savefig('example5notes.eps', format = 'eps')
pylab.show()
