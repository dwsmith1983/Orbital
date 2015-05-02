#!/usr/bin/env ipython
#  Solves example 6 of 230 Notes

import numpy as np
from scipy.integrate import odeint
import pylab
from scipy.optimize import brentq

me = 5.974e24  #  mass of the earth
mm = 7.348e22  #  mass of the moon
G = 6.67259e-20  #  gravitational parameter
re = 6378.0  #  radius of the earth in km
rm = 1737.0  #  radius of the moon in km
r12 = 384400.0  #  distance between the CoM of the earth and moon
d = 300  #  distance the spacecraft is above the Earth
pi1 = me / (me + mm)
pi2 = mm / (me + mm)
mue = 398600.0  #  gravitational parameter of earth km^3/sec^2
mum = G * mm  #  grav param of the moon
mu = mue + mum
omega = np.sqrt(mu / (r12 ** 3))

nu = -np.pi / 4  #  true anomaly  pick yourself

xl4 = r12 / 2 - pi2 * r12  #  x location of L4
yl4 = np.sqrt(3) / 2 * r12  #  y

print "The location of L4 is", xl4, yl4


#  Solve for Jacobi's constant
def f(C):
    return (omega ** 2 * (xl4 ** 2 + yl4 ** 2) + 2 * mue / r12 + 2 * mum / r12
            + 2 * C)


c = brentq(f, -5, 0)

print "Jacobi's constant is",c


x0 = (re + 200) * np.cos(nu) - pi2 * r12
#  x location of the satellite
y0 = (re + 200) * np.sin(nu)  #  y location

print "The satellite's initial position is", x0, y0

vbo = (np.sqrt(omega ** 2 * (x0 ** 2 + y0 ** 2) + 2 * mue /
               np.sqrt((x0 + pi2 * r12) ** 2 + y0 ** 2) + 2 * mum /
               np.sqrt((x0 - pi1 * r12) ** 2 + y0 ** 2) + 2 * -1.21))

print "Burnout velocity is", vbo

gamma = 0.4678 * np.pi / 180  #  flight path angle pick yourself

vx = vbo * (np.sin(gamma) * np.cos(nu) - np.cos(gamma) * np.sin(nu))
#  velocity of the bo in the x direction
vy = vbo * (np.sin(gamma) * np.sin(nu) + np.cos(gamma) * np.cos(nu))
#  velocity of the bo in the y direction

print "The satellite's initial velocity is", [vx, vy]

#  r0 = [x, y, 0]
#  v0 = [vx, vy, 0]
u0 = [x0, y0, 0, vx, vy, 0]


def deriv(u, dt):
    return [u[3],  #  dotu[0] = u[3]
            u[4],  #  dotu[1] = u[4]
            u[5],  #  dotu[2] = u[5]
            (2 * omega * u[4] + omega ** 2 * u[0] - mue * (u[0] + pi2 * r12) /
            np.sqrt(((u[0] + pi2 * r12) ** 2 + u[1] ** 2) ** 3) - mum *
            (u[0] - pi1 * r12) /
            np.sqrt(((u[0] - pi1 * r12) ** 2 + u[1] ** 2) ** 3)),
            #  dotu[3] = that
            (-2 * omega * u[3] + omega ** 2 * u[1] - mue * u[1] /
            np.sqrt(((u[0] + pi2 * r12) ** 2 + u[1] ** 2) ** 3) - mum * u[1] /
            np.sqrt(((u[0] - pi1 * r12) ** 2 + u[1] ** 2) ** 3)),
            #  dotu[4] = that
            0]  #  dotu[5] = 0

dt = np.linspace(0.0, 540000.0, 1000000.0)
u = odeint(deriv, u0, dt, rtol = 1e-13, atol = 1e-13)
x, y, z, vx, vy, vz = u.T

fig = pylab.figure()
ax = fig.add_subplot(111)
ax.plot(x, y, color = 'r')
#  adding the Lagrange point
L4 = pylab.Circle((xl4, yl4), radius = 500, color = '#000000')
ax.add_patch(L4)
#  adding the Earth
earth = pylab.Circle((-pi2 * r12, 0), radius = re, color = 'b')
ax.add_patch(earth)
#  adding the moon
moon = pylab.Circle((r12 * (1 - pi2), 0), radius = rm, color = 'gray')
ax.add_patch(moon)

pylab.xlim((-80000.0, 390000.0))
pylab.ylim((-80000.0, 390000.0))
pylab.savefig('example6notes.eps', format = 'eps')
pylab.show()

#  The code below finds the distance between path and l4
my_x, my_y, my_z = (xl4, yl4, 0.0)

delta_x = x - my_x
delta_y = y - my_y
delta_z = z - my_z
distance = np.array([np.sqrt(delta_x ** 2 + delta_y ** 2 + delta_z ** 2)])

minimum = np.amin(distance)
minimumtime = np.argmin(distance)


print "Distance from L4 at closest approach is", minimum
print "Location at closest approach is", [x[minimumtime], y[minimumtime]]

#  The code below finds the minimum speed at that location

vxmin = vx[minimumtime]  #  velocity in x at the min time
vymin = vy[minimumtime]  #  in y

vmin = np.sqrt(vxmin ** 2 + vymin **2)

print "The speed at the minimum distance is", vmin
print "The satellite's velocity at the minimum distance is", [vxmin, vymin]
