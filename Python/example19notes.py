#!/usr/bin/env ipython
#  This program solves example problem 19 of orbital 230 notes.

import numpy as np
import pylab
from scipy.optimize import fsolve
from numpy import linalg as LA
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
from matplotlib.ticker import MaxNLocator

Re = 1.496e8  #  semi-major axis of the Earth
Te = 365.25 * 24.0 * 3600.0  #  period of the Earth in sec
mus = 132712000000.0  #  grav param of the Sun


def f(a):
    return (2 * np.pi / np.sqrt(mus) * np.sqrt(a ** 3) - Te * 2.0 / 3.0)


a = fsolve(f, 100000000)

print 'The semi-major axis of the pre-flyby ellipse is ' + str(a)

e = Re / a - 1

print 'The eccentricity of the ellipse is ' + str(e)

vE = np.sqrt(mus / Re)

print 'The speed of the Earth is ' + str(vE)

rp = a * (1 - e)

print 'The radius at periapsis is ' + str(rp)

h = np.sqrt(2 * mus) * np.sqrt(Re * rp / (Re + rp))

print 'The angular momentum is ' + str(h)

vp = h / Re

print 'The velocity of the spacecraft of periapsis of the hyperbola is ' \
+ str(vp)

vinf = vE - vp

print 'The hyperbolic excess speed is ' + str(vinf)

alt = 500.0  #  the flyby distance
rph = 6378 + alt  #  radius at periapsis of the flyby hyperbola
mue = 398600.0  #  grav param of the Earth

eh = 1 + rph * vinf ** 2 / mue

print 'The eccentricity of the flyby hyperbola is ' + str(eh)

beta = np.arccos(1.0 / eh) * 180.0 / np.pi

print 'The angle at closest approach, beta, is ' + str(beta)

delta = 2.0 * beta

print 'The turn angle, delta, is ' + str(delta)

vpvec = np.array([0, -vp])

vEvec = np.array([0, -vE])

vinfoutvec = vinf * np.array([-np.sin(delta / 180.0 * np.pi),
                             np.cos(delta / 180.0 * np.pi)])

vhpostvec = np.array([vinfoutvec[0], vinfoutvec[1] + vpvec[1]])

vhpost = LA.norm(vhpostvec)

print 'The speed of the exiting velocity vector is ' + str(vhpost)


class FixedOrderFormatter(ScalarFormatter):
    """Formats axis ticks using scientific notation with a constant order of
    magnitude"""
    def __init__(self, order_of_mag = 0, useOffset = True, useMathText = False):
        self._order_of_mag = order_of_mag
        ScalarFormatter.__init__(self, useOffset = useOffset,
                                 useMathText = useMathText)
    def _set_orderOfMagnitude(self, range):
        """Over-riding this to avoid having orderOfMagnitude reset elsewhere"""
        self.orderOfMagnitude = self._order_of_mag


nu = np.linspace(-np.pi, np.pi, 500000)


def g(nu):
    return rph * (1.0 + eh) / (1.0 + eh * np.cos(nu - beta * np.pi / 180.0))


gt = g(nu)
ext = [np.argmin(gt), np.argmax(gt)]
gt[ext] = np.nan

fig = pylab.figure()
ax = fig.add_subplot(111)
ax.plot(gt * np.cos(nu), gt * np.sin(nu), 'r')

earth = pylab.Circle((0, 0), radius = 6378, color = 'b')
ax.add_patch(earth)

pylab.xlim((-50000, 20000))
pylab.ylim((-50000, 20000))

# Force the y-axis ticks to use 1e3 as a base exponent
ax.yaxis.set_major_formatter(FixedOrderFormatter(3))

ax.xaxis.set_major_formatter(FixedOrderFormatter(3))

pylab.gca().xaxis.set_major_locator(MaxNLocator(prune = 'lower'))

pylab.savefig('example19.eps', format = 'eps')

#pylab.show()



apost = 1.0 / (2.0 / Re - vhpost ** 2 / mus)

print 'The post hyperbolic flyby semi-major axis is ' + str(apost)

r0 = np.array([-Re, 0, 0])
v0 = np.array([vhpostvec[0], vhpostvec[1], 0])
h0vec = np.cross(r0, v0)
h0 = LA.norm(h0vec)

print 'The angular momentum of the new ellipse is ' + str(h0)

e2vec = np.cross(v0, h0vec) / mus - r0 / LA.norm(r0)

print 'The eccentricity vector of the new ellipse is ' + str(e2vec)

e2 = LA.norm(e2vec)

print 'The eccentricity is then ' + str(e2)

nu0 = np.arccos(np.dot(e2vec, np.array([1.0, 0, 0])) / e2)

nupost = (np.arccos(np.dot(r0, e2vec.ravel()) / (Re * e2)) * 180.0 / np.pi)

print 'The post flyby true anomaly is ' + str(nupost)
