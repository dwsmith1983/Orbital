#!/usr/bin/env ipython
#  This program solves example 9 from the 230 Orbital Notes,
#  Universal Variables
#  This program requires the use of stumpff.py

import numpy as np
from numpy import linalg as LA
import stumpff
from scipy.optimize import fsolve
import pylab
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
from matplotlib.ticker import MaxNLocator

class FixedOrderFormatter(ScalarFormatter):
    """Formats axis ticks using scientific notation with a constant order of
    magnitude"""
    def __init__(self, order_of_mag = 0, useOffset = True,
                 useMathText = False):
        self._order_of_mag = order_of_mag
        ScalarFormatter.__init__(self, useOffset = useOffset,
                                 useMathText = useMathText)
    def _set_orderOfMagnitude(self, range):
        """Over-riding this to avoid having orderOfMagnitude reset elsewhere"""
        self.orderOfMagnitude = self._order_of_mag

r0 = [7000.0, -12124.0, 0.0]  #  initial location
v0 = [2.6679, 4.6210, 0.0]    #  initial velocity
mue = 398600.0                #  grav param of Earth
deltat = 3600.0               #  change in time in sec
re = 6378.0                   #  radius of Earth in km

#  semi-major axis
a = 1 / (2.0 / LA.norm(r0) - LA.norm(v0) ** 2 / mue)

print 'The semi-major axis is ' + str(a)

sigma0 = np.dot(r0, v0) / np.sqrt(mue)

alpha = 1.0 / a

ecc = 1.0 / mue * np.cross(v0, np.cross(r0, v0)) - r0 / LA.norm(r0)

print 'The eccentric vector is ' + str(ecc)
print 'The eccentricity is ' + str(LA.norm(ecc))

#  Initial true anomaly in radians
nu0 = np.arccos(np.dot(r0, ecc) / (LA.norm(r0) * LA.norm(ecc)))

print 'The initial true anomaly is ' + str(nu0 / np.pi * 180) + ' degrees'

theta0 = np.arccos(np.dot([1.0, 0, 0], ecc) / LA.norm(ecc))


def f(x):
    return (sigma0 * x ** 2 * stumpff.C(alpha * x ** 2) +
            (1 - LA.norm(r0) * alpha) * x ** 3 * stumpff.S(alpha * x ** 2) +
            LA.norm(r0) * x - np.sqrt(mue) * deltat)


chi = fsolve(f, 200)

print 'The universal anomaly is ' + str(chi)

#  Final position vector
r = ((1 - chi ** 2 / LA.norm(r0) * stumpff.C(alpha * chi ** 2)) * r0 +
     (deltat - chi ** 3 / np.sqrt(mue) * stumpff.S(alpha * chi ** 2)) * v0)

print 'The final position vector is ' + str(r)

#  Final velocity vector
v = (chi * np.sqrt(mue) / (LA.norm(r) * LA.norm(r0)) *
     (alpha * chi ** 2 * stumpff.S(alpha * chi ** 2) - 1) * r0 +
     (1 - chi ** 2 / LA.norm(r) * stumpff.C(alpha * chi ** 2)) * v0)

print 'The final velocity vector is ' + str(v)

#  The final true anomaly is radians
nu = np.arccos(np.dot(r, ecc) / (LA.norm(r) * LA.norm(ecc)))

print 'The final true anomaly is ' + str(nu / np.pi * 180)


def g(theta):
    return ((-a * (LA.norm(ecc) ** 2 - 1)) /
            (1 + LA.norm(ecc) * np.cos(theta - theta0)))


theta = np.linspace(0, 2 * np.pi, 1000000)

fig = pylab.figure()
ax = fig.add_subplot(111, aspect = 'equal')
earth = pylab.Circle((0, 0), radius = re, color = 'b')
ax.plot(g(theta) * np.cos(theta), g(theta) * np.sin(theta), 'r')
ax.add_patch(earth)
ini = pylab.Circle((r0[0], r0[1]), radius = 150, color = '#000000')
fin = pylab.Circle((r[0], r[1]), radius = 150, color = '#000000')
ax.add_patch(ini)
ax.add_patch(fin)
pylab.xlim((-20000, 10000))
pylab.ylim((-20000, 10000))
ax.xaxis.set_major_formatter(FixedOrderFormatter(3))
ax.yaxis.set_major_formatter(FixedOrderFormatter(3))
ax.annotate('initial position', xy = (r0[0] * .3, r0[1] * 1.1))
ax.annotate('final position', xy = (r[0] * 1.7, r[1] * 1.1))
pylab.gca().xaxis.set_major_locator(MaxNLocator(prune = 'lower'))
pylab.savefig('example9notes.eps', format = 'eps')
pylab.show()
