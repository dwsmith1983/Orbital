#!/usr/bin/env ipython
#  This program solves part (a) and (b) of example 7 from the notes on
#  Bessel series approximation of the Kepler problem.

import numpy as np
from scipy.optimize import fsolve
import scipy.special as sp
import pylab

e = 0.0935  #  eccentricity of Mars
a = 227.9e6  #  semi-major axis of Mars
mus = 132712000000.0
time = 100  #  days

t = time * 24 * 3600  #  time in sec

energy = - mus / (2 * a)  #  energy of the ellipse

print 'The energy of the ellipse is ' + str(energy)

h = np.sqrt(-0.5 * mus ** 2 / energy * (1 - e ** 2))
#  angular momentum

print 'The angular momentum is ' + str(h)

Me = mus **2 / h ** 3 * (1 - e ** 2) ** 1.5 * t

print 'The mean anomaly is ' + str(Me)


def f(x):
    return x - e * np.sin(x) - Me


E = fsolve(f, .8)

print 'The eccentric anomaly is ' + str(E)

for m in range(1, 10, 1):
    E2 = (Me + sum(2.0 / n * sp.jv(n, n * e) * np.sin(n * Me)
             for n in range(1, m, 1)))
    print E2
