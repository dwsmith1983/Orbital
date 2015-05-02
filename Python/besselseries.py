#!/usr/bin/env ipython
#  This program plots the Bessel series

import numpy as np
import pylab as py
import scipy.special as sp

x = np.linspace(0, 15, 500000)


for v in range(0, 6):
    py.plot(x, sp.jv(v, x))


py.xlim((0, 15))
py.ylim((-0.5, 1.1))
py.legend(('$\mathcal{J}_0(x)$', '$\mathcal{J}_1(x)$', '$\mathcal{J}_2(x)$',
           '$\mathcal{J}_3(x)$', '$\mathcal{J}_4(x)$', '$\mathcal{J}_5(x)$'),
           loc = 0)
py.xlabel('$x$')
py.ylabel('$\mathcal{J}_n(x)$')
#py.title('Plots of the first six Bessel Functions')
py.grid(True)
#py.savefig('besseln0to6.eps', format = 'eps')


e = 0.99


def E(M):
    return (M + sum(2.0 / n * sp.jv(n, n * e) * np.sin(n * M)
                    for n in range(1, 3, 1)))

M = np.linspace(0, 2 * np.pi, 500000)

fig2 = py.figure()
ax2 = fig2.add_subplot(111, aspect = 'equal')
ax2.plot(E(M), M, 'b')


def E2(M):
    return (M + sum(2.0 / n * sp.jv(n, n * e) * np.sin(n * M)
                    for n in range(1, 11, 1)))


ax2.plot(E2(M), M, 'r')
py.xlim((0, 2 * np.pi))
py.ylim((0, 2 * np.pi))
py.xlabel('Eccentric anomaly, $E$')
py.ylabel('Mean anomaly, $M_e$')
py.legend(('$n = 3$', '$n = 10$'), loc = 0)
py.savefig('besselapproxn3and10.eps', format = 'eps')
py.show()
