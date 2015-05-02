#!/usr/bin/env ipython
#  This program plots the Stumpff functions C(z) and S(z)

import numpy as np
import pylab
from matplotlib.ticker import MaxNLocator


def C(z):
    if z > 0:
        return (1 - np.cos(z ** 0.5)) / z
    elif z < 0:
        return (np.cosh(np.sqrt(-z)) - 1) / -z
    return 0.5


def S(z):
    if z > 0:
        return (np.sqrt(z) - np.sin(z ** 0.5)) / np.sqrt(z) ** 3
    elif z < 0:
        return (np.sinh(np.sqrt(-z)) - np.sqrt(-z)) / np.sqrt(-z) ** 3
    return 1.0 / 6.0


if __name__ == '__main__':
    vC = np.vectorize(C)
    vS = np.vectorize(S)

    z = np.linspace(-50.0, 500.0, 100000.0)
    y = vC(z)
    y2 = vS(z)


    fig = pylab.figure()
    ax = fig.add_subplot(111)
    ax.plot(z, y, 'r')
    ax.plot(z, y2, 'b')
    pylab.legend(('$C(z)$', '$S(z)$'), loc = 0)
    pylab.xlim((-50, 0))
    pylab.ylim((0, 12))
    pylab.xlabel('$z$')
    pylab.gca().xaxis.set_major_locator(MaxNLocator(prune = 'lower'))
    pylab.savefig('stumpffneg50to0.eps', format = 'eps')


    fig2 = pylab.figure()
    ax2 = fig2.add_subplot(111)
    ax2.plot(z, y, 'r')
    ax2.plot(z, y2, 'b')
    pylab.legend(('$C(z)$', '$S(z)$'), loc = 1)
    pylab.xlim((0, 30))
    pylab.ylim((0, 0.5))
    pylab.xlabel('$z$')
    pylab.gca().xaxis.set_major_locator(MaxNLocator(prune = 'lower'))
    pylab.savefig('stumpff0to30.eps', format = 'eps')


    fig3 = pylab.figure()
    ax3 = fig3.add_subplot(111)
    ax3.plot(z, y, 'r')
    ax3.plot(z, y2, 'b')
    pylab.legend(('$C(z)$', '$S(z)$'), loc = 0)
    pylab.xlim((0, 500))
    pylab.ylim((0, 0.05))
    pylab.xlabel('$z$')
    pylab.gca().xaxis.set_major_locator(MaxNLocator(prune = 'lower'))
    pylab.savefig('stumpff0to500.eps', format = 'eps')
    pylab.show()
