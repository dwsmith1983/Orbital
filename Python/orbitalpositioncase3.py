#!/usr/env/bin ipython
#  Plots the Eccentric anomaly as a function of the true anomaly

import pylab
import numpy as np


def E(nup):
    return .5 * np.tan(nu / 2) + 1 / 6 * np.tan(nu / 2) ** 3


nu = np.linspace(0.0, 2.0 * np.pi, 50000.0)

E2 = E(nu)
ext = [np.argmin(E2), np.argmax(E2)]
E2[ext] = np.nan


fig = pylab.figure()
ax = fig.add_subplot(111)

ax.plot(nu, E2)

pylab.legend()
pylab.xlim((0, np.pi))
pylab.ylim((0, np.pi))
#pylab.savefig('eccentruepara.eps', format = 'eps', dpi = 1000)
pylab.show()
