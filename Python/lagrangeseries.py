#!/usr/bin/env ipython
#  plots the Lagrange series of the mean anomaly

import numpy as np
import pylab as py
from scipy.misc import factorial as fact

e = 0.65


def E(M):
    return (M + sum((1.0 / 2.0 ** (n - 1) * sum((-1) ** k /
                                                (fact(n - k) * fact(k)) *
                                                (n - 2 * k) ** (n - 1) *
                                                np.sin((n - 2 * k) * M)
                                                for k in range(0, n / 2, 1))) *
                                                e ** n
                                                for n in range(1, 4, 1)))


M = np.linspace(0, 2 * np.pi, 50000.0)

fig = py.figure()
ax = fig.add_subplot(111, aspect = 'equal')
ax.plot(E(M), M, 'b')
#py.savefig('lagrangeeccen65n3.eps', format = 'eps')
#py.show()


def E2(M):
    return (M + sum((1.0 / 2.0 ** (n - 1) * sum((-1) ** k /
                                                (fact(n - k) * fact(k)) *
                                                (n - 2 * k) ** (n - 1) *
                                                np.sin((n - 2 * k) * M)
                                                for k in range(0, n / 2, 1))) *
                                                e ** n
                                                for n in range(1, 11, 1)))


M = np.linspace(0, 2 * np.pi, 50000.0)

#fig2 = py.figure()
#ax2 = fig2.add_subplot(111)
ax.plot(E2(M), M, 'r')
ax.legend(('n = 3', 'n = 10'), loc = 0)
py.axis([0, 2 * np.pi, 0, 2 * np.pi])
py.savefig('lagrangeeccen65n3n10.eps', format = 'eps')
py.show()


e = 0.9


def E3(M):
    return (M + sum((1.0 / 2.0 ** (n - 1) * sum((-1) ** k /
                                                (fact(n - k) * fact(k)) *
                                                (n - 2 * k) ** (n - 1) *
                                                np.sin((n - 2 * k) * M)
                                                for k in range(0, n / 2, 1))) *
                                                e ** n
                                                for n in range(1, 4, 1)))


M = np.linspace(0, 2 * np.pi, 50000.0)

fig2 = py.figure()
ax2 = fig2.add_subplot(111, aspect = 'equal')
ax2.plot(E3(M), M, 'b')
#py.savefig('lagrangeeccen65n3.eps', format = 'eps')
#py.show()


def E4(M):
    return (M + sum((1.0 / 2.0 ** (n - 1) * sum((-1) ** k /
                                                (fact(n - k) * fact(k)) *
                                                (n - 2 * k) ** (n - 1) *
                                                np.sin((n - 2 * k) * M)
                                                for k in range(0, n / 2, 1))) *
                                                e ** n
                                                for n in range(1, 11, 1)))


M = np.linspace(0, 2 * np.pi, 50000.0)

#fig2 = py.figure()
#ax2 = fig2.add_subplot(111)
ax2.plot(E4(M), M, 'r')
py.xlim((0, 2 * np.pi))
py.ylim((0, 2 * np.pi))
ax2.legend(('n = 3', 'n = 10'), loc = 0)
py.axis([0, 2 * np.pi, 0, 2 * np.pi])
py.savefig('lagrangeeccen90n3n10.eps', format = 'eps')
py.show()
