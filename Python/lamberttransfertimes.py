#!/usr/bin/env ipython
#  This program plots the different transfer times vs a.

import numpy as np
import pylab

r1 = 1  #  AU Earth
r2 = 1.524  #  AU Mars
deltanu = 75 * np.pi / 180  #  angle in radians
mu = 38.26836055644185

c = np.sqrt(r1 ** 2 + r2 ** 2 - 2 * r1 * r2 * np.cos(deltanu))

print 'The chord of the space triangle is ' + str(c)

s = (r1 + r2 + c) / 2

print 'The semi-perimeter of the space triangle is ' + str(s)

am = s / 2

print 'The minimum semi-major axis is ' + str(am)

tm = (np.sqrt(s ** 3 / 8) * (np.pi - 2 * np.arcsin(np.sqrt(1 - c / s)) +
                            np.sin(2 * np.arcsin(np.sqrt(1 - c / s)))))

days = tm * 365.25 / (2 * np.pi)

print 'The minimum transfer time is ' + str(days) + ' days.'


def g(a):
    alphag = 2* np.pi - 2 * np.arcsin(np.sqrt(s / (2 * a)))
    betag = 2 * np.arcsin(np.sqrt((s - c) / (2 * a)))
    return (np.sqrt(a ** 3 / mu)
            * (alphag - betag - (np.sin(alphag) - np.sin(betag))))


def f(a):
    alpha = 2 * np.arcsin(np.sqrt(s / (2 * a)))
    beta = 2 * np.arcsin(np.sqrt((s - c) / (2 * a)))
    return (np.sqrt(a **3 / mu) * (alpha - beta - (np.sin(alpha)
                                                      - np.sin(beta))))


a = np.linspace(am, 200, 500000)

fig = pylab.figure()
ax = fig.add_subplot(111)
ax.plot(a, f(a), color = '#000000')
ax.plot(a, g(a), color = '#000000')
ax.plot([am, am], [f(am), g(am)], color = '#000000')
ax.plot([am, am], [0, 0.5039140423825622], 'b--')
ax.plot([0, am], [0.5039140423825622, 0.5039140423825622], 'b--')
ax.plot([0, 200], [0.197, 0.197], 'b--')
pylab.xlim((0.9, 2))
pylab.ylim((0, 2))

pylab.xlabel('Semi-major Axis $a$ in AU')
pylab.ylabel('Time of Flight in Years')
pylab.text(1, 1.7, '$r_1 = 1.0$ AU\n' +\
           '$r_2 = 1.524$ AU\n' +\
           '$\\Delta \\nu = 75^{\\circ}$',
           fontsize = 11, color = 'r')
pylab.text(1.6, 0.3, '$\\alpha = \\alpha_0$', fontsize = 14,
           color = 'r')
pylab.text(0.5, 0.9, '$\\alpha = 2\\pi - \\alpha_0$',
           transform = ax.transAxes, rotation = 40, fontsize = 14,
           color = 'r')
pylab.savefig('lamberttransfertimes.eps', format = 'eps')
pylab.show()
