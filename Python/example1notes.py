#!/usr/bin/env ipython
#  Example 1 of the 230 Orbital Notes--Elliptical

import numpy as np
import pylab

e = np.arange(0, 1, 0.08).reshape(-1, 1)
nu = np.linspace(0, 2 * np.pi, 5000).reshape(1, -1)
vvals = np.sqrt((e ** 2) * np.ones(nu.shape) + 2 * e * np.cos(nu) + 1)

for v, _e in zip(vvals, e.ravel()):
    pylab.plot(nu.ravel(), v, label = str(_e))

pylab.legend()
pylab.xlim(0, 8)
pylab.savefig('example1notes.eps', format = 'eps', dpi = 1000)
pylab.show()
