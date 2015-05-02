#!/usr/bin/env ipython
#  This program solves problem 18 in 230 Orbital notes

import numpy as np

mus = 1.327e11  #  grav param of the Sun
mum = 42830.0  #  grav param of Mars
Re = 1.496e8  #  semi-major axis Earth
Rm = 2.279e8  #  semi-major axis Mars
rm = 3396.0  #  the radius of Mars
T = 7.0 * 3600.0  #  period in sec

vinf = np.sqrt(mus / Rm) * (1 - np.sqrt(2 * Re / (Re + Rm)))

print 'The hyperbolic excess speed on arrival to Mars is ' + str(vinf)

a = (T * np.sqrt(mum) / (2 * np.pi)) ** (2.0 / 3.0)

print 'The semi-major axis of the capture orbit is ' + str(a)

e = 2 * mum / (a * vinf ** 2) - 1

print 'The eccentricity of the capture orbit is ' + str(e)

deltaV = vinf * np.sqrt((1 - e) / 2)

print 'The minimum delta-v requirement is ' + str(deltaV)

rp = ((2 * mum) * (1 - e)) / (vinf ** 2 * (1 + e))

print 'The radius at periapsis is ' + str(rp)

delta = rp * np.sqrt(2 / (1 - e))

print 'The aiming radius is ' + str(delta)

beta = np.arccos(1 / (1 + rp * vinf ** 2 / (mum))) * 180 / np.pi

print 'The angle between periapsis and Mars velocity vector is ' + str(beta)
