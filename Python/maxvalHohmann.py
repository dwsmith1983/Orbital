#!/usr/bin/env ipython
#  This programs plots and finds the max value of the function
#  \Delta\bar{v}_H

import pylab
import numpy as np
from scipy.optimize import fmin


def f(R):
    return (1.0 / np.sqrt(R) - np.sqrt(2) * (1 - R) /
            np.sqrt(R * (1 + R)) - 1)

R = np.linspace(1, 20, 500000)

max_r = fmin(lambda r: 1.0 / f(r), 20)

fig = pylab.figure()
ax = fig.add_subplot(111)
ax.plot(R, f(R), color = 'r')
pylab.xlabel('$R_1 = \\frac{r_C}{r_A}$')
pylab.savefig('maxvalHohmann.eps', format = 'eps')


print 'The maximum value of f(R) is ' + str(max_r)


def g(x, r):
    return (np.sqrt(2.0 * (x + r) / (r * x)) - (1 + np.sqrt(x)) / np.sqrt(x) -
            np.sqrt(2.0 / (r * (1 + r))) * (1 - r) -
            (1.0 / np.sqrt(x) - np.sqrt(2.0) * (1 - x) / np.sqrt(x * (1 + x))
             - 1))


x, r = np.mgrid[1:50:200j, 1:100:200j]
Z = g(x, r)

fig2 = pylab.figure()
ax2 = fig2.add_subplot(111)
ax2.contour(x, r, Z, colors = 'r', levels = [0])
ax2.plot([11.94, 11.94], [0, 1000], 'k--')
pylab.xlabel('$R_1 = \\frac{r_C}{r_A}$')
pylab.ylabel('$R_2 = \\frac{r_B}{r_A}$')
pylab.xlim((0, 35))
pylab.ylim((0, 100))
pylab.text(20, 17, '$r_B = r_C$', fontsize = 11, color = 'k')
pylab.text(20, 60, 'Bi-elliptic', fontsize = 11, color = 'k')
pylab.text(5, 60, 'Hohmann', fontsize = 11, color = 'k')
pylab.savefig('biellpvshoh.eps', format = 'eps')

pylab.show()
