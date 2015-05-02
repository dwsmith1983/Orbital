#!/usr/env/bin ipython
#  Plots the Eccentric anomaly as a function of the true anomaly

import pylab
import numpy as np

e = np.arange(0.0, 1.0, 0.15).reshape(-1, 1)
E = np.linspace(0, 2 * np.pi, 50000).reshape(1, -1)
Mevals = E - e * np.sin(E)

fig = pylab.figure()
ax = fig.add_subplot(111)

for Me, _e in zip(Mevals, e.ravel()):
    ax.plot(E.ravel(), Me, label = str(_e))

pylab.legend(loc = 'upper left')
pylab.xlim((0.0, 2 * np.pi))
pylab.ylim((0, 2 * np.pi))
#pylab.savefig('eccenmeananomfunc.eps', format = 'eps', dpi = 1000)
pylab.show()

nu = np.linspace(0.001, 2 * np.pi - 0.001, 50000)
x = ((1 - e) / (1 + e)) ** 0.5 * np.tan(nu / 2)
M2evals = (2 * np.arctan2(1, 1 / x) -
           e * (1 - e ** 2) ** 0.5 * np.sin(nu) / (1 + e * np.cos(nu)))

fig2 = pylab.figure()
ax2 = fig2.add_subplot(111)

for Me2, _e in zip(M2evals, e.ravel()):
    ax2.plot(nu.ravel(), Me2, label = str(_e))

pylab.legend(loc = 'upper left')
pylab.xlim((0, 2 * np.pi))
pylab.ylim((0, 2 * np.pi))
#pylab.savefig('eccentrueanomfunc.eps', format = 'eps', dpi = 1000)
pylab.show()
