#!/usr/bin/env ipython
#  This program plots two elliptical orbits in different planes

import numpy as np
from scipy.integrate import odeint
from mpl_toolkits.mplot3d import Axes3D
import pylab

mu = 1.0  #  grav param of the Sun in AU

u0 = [1.0, 0, 0, .25, 1.0, 0]
w0 = [4.25, 0, 0, .5, 2.25, 1.0]


def deriv(u, t):
    n = -mu / np.sqrt(u[0]**2 + u[1]**2 + u[2]**2)
    return [u[3],      # u[0]' = u[3]
            u[4],      # u[1]' = u[4]
            u[5],      # u[2]' = u[5]
            u[0] * n,  # u[3]' = u[0] * n
            u[1] * n,  # u[4]' = u[1] * n
            u[2] * n]  # u[5]' = u[2] * n


#def deriv2(w, t):
#    n = -mu / np.sqrt(w[0]**2 + w[1]**2 + w[2]**2)
#    return [w[3],      # u[0]' = u[3]
#            w[4],      # u[1]' = u[4]
#            w[5],      # u[2]' = u[5]
#            w[0] * n,  # u[3]' = u[0] * n
#            w[1] * n,  # u[4]' = u[1] * n
#            w[2] * n]  # u[5]' = u[2] * n


dt = np.linspace(0.0, 2 * np.pi + .25, 200000.0)
#deltt = np.linspace(0.0, 51712.2, 200000.0)


u = odeint(deriv, u0, dt)
x, y, z, x1, y1, z1 = u.T
#w = odeint(deriv2, w0, dt, atol = 1e-14, rtol = 1e-14)
#x2, y2, z2, x3, y3, z3 = u.T

fig = pylab.figure()
ax = fig.add_subplot(111, projection = '3d', aspect = 'equal')
ax.plot(x, y, z, 'r')
#ax.plot(x2, y2, z2, 'k')

pylab.show()
