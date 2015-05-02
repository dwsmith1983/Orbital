#!/usr/env/bin ipython
#  Plots the Eccentric anomaly as a function of the true anomaly

import pylab
import numpy as np

e = np.arange(1.1, 5.0, 0.4).reshape(-1, 1)
nu = np.linspace(0.01, 2 * np.pi, 50000).reshape(1, -1)
Mhvals = (e * np.sin(nu) * (e ** 2 - 1) ** 0.5 / (1 + e * np.cos(nu))
           - np.log(((e + 1) ** 0.5 + (e - 1) ** 0.5 * np.tan(nu / 2)) /
                    ((e + 1) ** 0.5 - (e - 1) ** 0.5 * np.tan(nu / 2))))


fig = pylab.figure()
ax = fig.add_subplot(111)

for Mh, _e in zip(Mhvals, e.ravel()):
    ax.plot(nu.ravel(), Mh, label = str(_e))

pylab.legend(loc = 'upper left')
pylab.xlim((0, np.pi))
pylab.ylim((0, 10))
#pylab.savefig('eccentruehyp.eps', format = 'eps', dpi = 1000)
pylab.show()

F = np.linspace(0, 2 * np.pi, 50000)
Mh2vals = e * np.sinh(F) - F

fig2 = pylab.figure()
ax2 = fig2.add_subplot(111)

for Mh2, _e in zip(Mh2vals, e.ravel()):
    ax2.plot(F.ravel(), Mh2, label = str(_e))

pylab.legend(loc = 'upper left')
pylab.xlim((0, 2 * np.pi))
pylab.ylim((0, 800))
#pylab.savefig('eccenmeanhyp.eps', format = 'eps', dpi = 1000)
pylab.show()
