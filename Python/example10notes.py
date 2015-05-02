#!/usr/bin/env ipython
#  This program solves example 10 of 230 Orbital notes.

import numpy as np
import pylab
from matplotlib.ticker import MaxNLocator
import numpy.linalg as LA

r1 = 1.0  #  Earth AU
r2 = 5.2  #  Jupiter AU
mu = 1.0

def a(nu):
    return ((r1 + r2 + np.sqrt(r1 ** 2 + r2 ** 2 -
                               2 * r1 * r2 * np.cos(nu))) / 4)


nu = np.linspace(0, np.pi, 500000)

fig = pylab.figure()
ax = fig.add_subplot(111)
ax.plot(nu, a(nu), 'r')
pylab.gca().xaxis.set_major_locator(MaxNLocator(prune = 'lower'))
pylab.xlim((0, np.pi))
pylab.ylim((2.55, np.pi))

pylab.xlabel('$\\Delta\\nu$ in radians')
pylab.ylabel('$a_m$ in AU')
pylab.savefig('amfunctionofdeltnu.eps', format = 'eps')


nu2 = 150.0 / 180.0 * np.pi  #  change in true anomaly
a = 5  #  semi-major axis
am = (r1 + r2 + np.sqrt(r1 ** 2 + r2 ** 2 - 2 * r1 * r2 * np.cos(nu2))) / 4

print 'The minimum semi-major axis for delta nu = 150 is ' + str(am)

c = np.sqrt(r1 ** 2 + r2 ** 2 - 2 * r1 * r2 * np.cos(nu2))

print 'The chord of the space triangle is ' + str(c)

s = (r1 + r2 + c) / 2

print 'The semi-perimeter of the space triangle is ' + str(s)

betam = 2 * np.arcsin(np.sqrt(1 - c / s))
#  bm pos since deltnu < pi

tm = np.sqrt(s ** 3 / 8) * (np.pi - betam + np.sin(betam))

print 'The minimum energy time of flight is ' + str(tm) + ' time units.'

DU = 1.496e8  #  distance units
TU = DU ** 1.5 / np.sqrt(132712000000.0) / (24 * 3600)
#print TU
#  time units

print 'The minimum energy time of flight in days is ' + str(TU * tm)

#  short time of flight
alpha1 = 2 * np.arcsin(np.sqrt(s / (2 * a)))
beta1 = 2 * np.arcsin(np.sqrt((s - c) / (2 * a)))
t1 = a ** 1.5 * (alpha1 - beta1 - np.sin(alpha1) + np.sin(beta1))

print 'The short time of flight in days is ' + str(TU * t1)

#  long time of flight
alpha2 = 2 * np.pi - 2 * np.arcsin(np.sqrt(s / (2 * a)))
beta2 = beta1
t2 = a ** 1.5 * (alpha2 - beta2 - np.sin(alpha2) + np.sin(beta2))

print 'The long time of flight in days is ' + str(TU * t2)

r1v = np.array([1, 0])  #  r1 as a vector
r2v = 5.2 * np.array([np.cos(nu2), np.sin(nu2)])  #  r2 as a vector

u1 = r1v  #  unit vector of r1v
u2 = r2v / LA.norm(r2v)  #  unit vector of r2v
uc = (r2v - r1v) / c  #  unit vector along the chord

print 'The unit vector to Jupiter is ' + str(u2)
print 'The unit vector of the chord is ' + str(uc)

#  vector coefficients for the short flight:
As = np.sqrt(mu / (4.0 * a)) * 1.0 / np.tan(alpha1 / 2)
Bs = np.sqrt(mu / (4.0 * a)) * 1.0 / np.tan(beta1 / 2)

#print str(As)
#print str(Bs)

Al = np.sqrt(mu / (4.0 * a)) * 1.0 / np.tan(alpha2 / 2)
Bl = np.sqrt(mu / (4.0 * a)) * 1.0 / np.tan(beta2 / 2)

#print str(Al)
#print str(Bl)

v1s = (Bs + As) * uc + (Bs - As) * u1
v2s = (Bs + As) * uc + (As - Bs) * u2

print 'The velocity vector v1 for the short time is ' + str(v1s)
print 'The velocity vector v2 for the short time is ' + str(v2s)

v1l = (Bl + Al) * uc + (Bl - Al) * u1
v2l = (Bl + Al) * uc + (Al - Bl) * u2

print 'The velocity vector v1 for the short time is ' + str(v1l)
print 'The velocity vector v2 for the short time is ' + str(v2l)

eccs = (LA.norm(v1s) ** 2 - 1) * r1v - np.dot(r1v, v1s) * v1s

print 'The eccentric vector for the short time of flight is ' + str(eccs)

print 'The eccentricity for the short time of flight is ' + str(LA.norm(eccs))

ps = ((4 * a * (s - LA.norm(r1v)) * (s - LA.norm(r2v))) / c ** 2 *
      np.sin((alpha1 + beta1) / 2) ** 2)
#  semi-latus rectum short time of flight

print 'The semi-latus rectum for the short time of flight is ' + str(ps)

nuoffs = np.arccos(np.dot(r1v, eccs) / (LA.norm(r1v) * LA.norm(eccs)))
#  true anomaly offset for the eccentric vector

vectors = np.array([ [0, 0, r1v[0], r1v[1]], [0, 0, r2v[0], r2v[1]] ])
X, Y, U, V = zip(*vectors)


def t(theta):
    return ps / (1 + LA.norm(eccs) * np.cos(theta - nuoffs))


theta = np.linspace(0, 2 * np.pi, 500000)

fig2 = pylab.figure()
ax2 = fig2.add_subplot(111, aspect = 'equal')
ax2.plot(t(theta) * np.cos(theta), t(theta) * np.sin(theta), 'r')
ax2.plot(5.2 * np.cos(theta), 5.2 * np.sin(theta), 'b--')
ax2.plot(np.cos(theta), np.sin(theta), 'b--')
ax2.quiver(X, Y, U, V, angles = 'xy', scale_units = 'xy', scale = 1,
           width = 0.005)

pylab.xlim((-10, 6))
pylab.ylim((-8, 8))
pylab.savefig('example10short.eps', format = 'eps')

#  long time of flight

eccl = (LA.norm(v1l) ** 2 - 1) * r1v - np.dot(r1v, v1l) * v1l

print 'The eccentric vector for the long time of flight is ' + str(eccl)


print 'The eccentricity for the long time of flight is ' + str(LA.norm(eccl))

pl = ((4 * a * (s - LA.norm(r1v)) * (s - LA.norm(r2v))) / c ** 2 *
      np.sin((alpha2 + beta2) / 2) ** 2)
#  semi-latus rectum long time of flight

print 'The semi-latus rectum for the short time of flight is ' + str(pl)

nuoffl = np.arccos(np.dot(r1v, eccl) / (LA.norm(r1v) * LA.norm(eccl)))
#  true anomaly offset for the eccentric vector


def y(theta):
    return pl / (1 + LA.norm(eccl) * np.cos(theta + nuoffl))


fig3 = pylab.figure()
ax3 = fig3.add_subplot(111, aspect = 'equal')
ax3.plot(y(theta) * np.cos(theta), y(theta) * np.sin(theta), 'r')
ax3.plot(5.2 * np.cos(theta), 5.2 * np.sin(theta), 'b--')
ax3.plot(np.cos(theta), np.sin(theta), 'b--')
ax3.quiver(X, Y, U, V, angles = 'xy', scale_units = 'xy', scale = 1,
           width = 0.005)

pylab.xlim((-8, 7))
pylab.ylim((-6, 9))
pylab.savefig('example10long.eps', format = 'eps')

pylab.show()
