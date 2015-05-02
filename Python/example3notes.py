#!/usr/bin/env ipython
#  Example 3 of 230 Notes hyperbolic trajectory

import numpy as np
import pylab
from scipy.optimize import fsolve
from scipy.integrate import odeint
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


vin = 2.23  #  meteors incoming velocity
mue = 398600.0  #  grav param of Earth
r = 402000.0  #  the distance from Earth at t = 0
nu = 5.0 * np.pi / 6.0  #  true anomaly
re = 6378.0  #  the radius of the Earth

energy = vin ** 2.0 / 2.0 - mue / r

print "The energy of the incoming meteor is", energy

def f(e):
    return (mue ** 2 / (2.0 * energy) * (e ** 2.0 - 1.0) - r * mue *
            (1.0 - e * np.sqrt(3.0) / 2.0))

e = fsolve(f, 1)

print "The eccentricity of the trajectory is", e

h = np.sqrt(mue ** 2.0 * (e ** 2.0 - 1.0) / (2.0 * energy))
#  angular momentum

print "The angular momentum is", h

rp = h ** 2.0 / (mue * (1.0 + e))  #  the radius at periapsis

print "The radius at periapsis is", rp

alt = rp - re  #  the altitude at periapsis

print "The altitude at periapsis is", alt

v = h / rp  #  speed at closest approach

print "The speed at closest approach is", v

p = rp * (1.0 + e)  #  semi-latus rectum


def rr(theta):
    return p / (1.0 + e * np.cos(theta))


theta = np.linspace(0, 2 * np.pi, 25000 * np.pi)

#  removing the asymptotes
rt = rr(theta)
ext = [np.argmin(rt), np.argmax(rt)]
rt[ext] = np.nan

fig = pylab.figure()
ax = fig.add_subplot(111, aspect = 'equal')
ax.plot(rt * np.cos(theta), rt * np.sin(theta), 'r')
pylab.xlim((-40000.0, 40000.0))
pylab.ylim((-40000.0, 40000.0))
#  adding the earth
earth = pylab.Circle((0, 0), radius = re, color = 'b')
ax.add_patch(earth)
ax.xaxis.set_major_formatter(FixedOrderFormatter(3))
ax.yaxis.set_major_formatter(FixedOrderFormatter(3))
pylab.savefig('example3notes.eps', format = 'eps')
pylab.show()

vinf = mue / h * np.sqrt(e ** 2.0 - 1)  #  hyperbolic excess speed

print "The hyperbolic excess speed is", vinf

rx = r * np.cos(nu)  #  the x initial position
ry = r * np.sin(nu)  #  the y
vx = mue / h * np.sin(nu)  #  the vx initial speed
vy = -mue / h * (e + np.cos(nu))  #  the vy

print "The initial position of the meteor is", [rx, ry]
print "The initial velocity is", [vx, vy]

u0 = [rx, ry, 0.0, vx, vy, 0.0]

print "The initial vector is", u0


def deriv(u, dt):
    return [u[3],
            u[4],
            u[5],
            -mue * u[0] / (u[0] ** 2 + u[1] ** 2 + u[2] ** 2) ** 1.5,
            -mue * u[1] / (u[0] ** 2 + u[1] ** 2 + u[2] ** 2) ** 1.5,
            -mue * u[2] / (u[0] ** 2 + u[1] ** 2 + u[2] ** 2) ** 1.5]


dt = np.linspace(0.0, 500000.0, 1000000.0)
u = odeint(deriv, u0, dt, atol = 1e-13, rtol = 1e-13)
x, y, z, vx, vy, vz = u.T

#  find the closet approach
myx, myy = (0, 0)

deltax = x - myx
deltay = y - myy

distance = np.array([np.sqrt(deltax ** 2 + deltay ** 2)])

minimum = np.amin(distance)
print "The time to minimum distance is", np.argmin(distance)
print "The distance at periapsis is", minimum

fig2 = pylab.figure()
ax2 = fig2.add_subplot(111, aspect = 'equal')
ax2.plot(x, y, 'r')
pylab.xlim((-400000.0, 30000.0))
pylab.ylim((-205000.0, 205000.0))
#  adding the earth
earth2 = pylab.Circle((0,0), radius = re, color = 'b')
ax2.add_patch(earth2)
ax2.xaxis.set_major_formatter(FixedOrderFormatter(3))
ax2.yaxis.set_major_formatter(FixedOrderFormatter(3))
pylab.savefig('example3notes2.eps', format = 'eps')
pylab.show()
