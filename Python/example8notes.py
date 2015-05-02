#!/usr/bin/env ipython
#  This program solves the hyperbolic Kepler equation from example 7
#  of the notes.

import numpy as np
import pylab
from scipy.integrate import odeint
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter

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


a = -35000  #  the semi-major axis of the hyperbolic trajectory
e = 1.2  #  the eccentricity
nu0 = -20 * np.pi / 180  #  initial anomaly
nuf = 103 * np.pi / 180  #  final anomaly

F0 = 2 * np.arctanh(np.sqrt((e - 1) / (e + 1)) * np.tan(nu0 / 2))

print 'The initial eccentric anomaly is ' + str(F0)

Ff = 2 * np.arctanh(np.sqrt((e - 1) / (e + 1)) * np.tan(nuf / 2))

print 'The final eccentric anomaly is ' + str(Ff)

M0 = e * np.sinh(F0) - F0

print 'The initial mean anomaly is ' + str(M0)

Mf = e * np.sinh(Ff) - Ff

print 'The initial mean anomaly is ' + str(Mf)

n = np.sqrt(398600.0 / -a ** 3)  #  the mean motion

print 'The mean motion is ' + str(n)

deltat = (Mf - M0) / n  #  change of time in secs

hours = deltat / 3600.0

print 'The change in time from nu_0 to nu_f is ' + str(deltat) + ' sec'

print 'Or ' + str(hours) + ' hrs'

energy = -398600.0 / (2.0 * a)
h = np.sqrt(a * 398600.0 * (1 - e ** 2))


def r(nu):
    return h ** 2.0 / (398600.0 * (1.0 + e * np.cos(nu)))


nu = np.linspace(0.0, 2.0 * np.pi, 500000.0)

rt = r(nu)
ext = [np.argmin(rt), np.argmax(rt)]
rt[ext] = np.nan

fig = pylab.figure()
ax = fig.add_subplot(111, aspect = 'equal')
earth = pylab.Circle((0, 0), radius = 6378, color = 'blue')
ax.add_patch(earth)
ax.plot(rt * np.cos(nu), rt * np.sin(nu), 'r')
pylab.axis([-70000, 10000, -40000, 40000])
pylab.gca().xaxis.set_major_locator(MaxNLocator(prune = 'lower'))
ax.xaxis.set_major_formatter(FixedOrderFormatter(3))
ax.yaxis.set_major_formatter(FixedOrderFormatter(3))
pylab.savefig('example8notes.eps', format = 'eps')
pylab.show()


r = h ** 2 / (398600 * (1 + e * np.cos(nu0)))

rx = r * np.cos(nu0)
ry = r * np.sin(nu0)
r0 = [rx, ry]
vx = -398600.0 / h * np.sin(nu0)
vy = 398600.0 / h * (e + np.cos(nu0))
v0 = [vx, vy]

rf = h ** 2 / (398600 * (1 + e * np.cos(nuf)))

rxf = rf * np.cos(nuf)
ryf = rf * np.sin(nuf)


print 'The initial position of the spacecraft is ' + str(r0)
print 'The initial velocity of the spacecraft is ' + str(v0)

u0 = [rx, ry, 0, vx, vy, 0]

mue = 398600


def deriv(u, dt):
    return [u[3],
            u[4],
            u[5],
            -mue * u[0] / (u[0] ** 2 + u[1] ** 2 + u[2] ** 2) ** 1.5,
            -mue * u[1] / (u[0] ** 2 + u[1] ** 2 + u[2] ** 2) ** 1.5,
            -mue * u[2] / (u[0] ** 2 + u[1] ** 2 + u[2] ** 2) ** 1.5]


dt = np.linspace(0.0, deltat, 10.0)
u = odeint(deriv, u0, dt, atol = 1e-13, rtol = 1e-13)
x, y, z, vx, vy, vz = u.T

fig2 = pylab.figure()
ax2 = fig2.add_subplot(111, aspect = 'equal')
ax2.plot(x, y, 'ro')
pylab.xlim((-25000.0, 10000.0))
pylab.ylim((-10000.0, 25000.0))
ax2.set_xticks([5000.0])
#  adding the earth
earth2 = pylab.Circle((0,0), radius = 6378, color = 'b')
ax2.add_patch(earth2)
pylab.gca().xaxis.set_major_locator(MaxNLocator(prune = 'lower'))
ax2.xaxis.set_major_formatter(FixedOrderFormatter(3))
ax2.yaxis.set_major_formatter(FixedOrderFormatter(3))
pylab.savefig('example8notes2.eps', format = 'eps')
pylab.show()
