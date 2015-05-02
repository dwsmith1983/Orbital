#!/usr/bin/env ipython
#  This program solves example 13 for a Hohmann transfer and then
#  a phasing maneuver

import numpy as np

r1 = 8000.0
r2 = 14000.0
mue = 398600.0
re = 6378.0

v1 = np.sqrt(mue / r1)  #  speed of sat 1

print 'The speed of the first satellite is ' + str(v1)

v2 = np.sqrt(mue / r2)  #  speed of sat 2

print 'The speed of the second satellite is ' + str(v2)

a3 = .5 * (r1 + r2)  #  semi-major axis of the Hohmann ellipse

print 'The semi-major axis of the Hohmann ellipse is ' + str(a3)

tofsec = np.pi / np.sqrt(mue) * np.sqrt(a3 ** 3)

print 'The time of flight is ' + str(tofsec) + ' secs'

h3 = np.sqrt(2.0 * mue) * np.sqrt(r1 * r2 / (r1 + r2))
#  angular momentum of Hohmann

vp3 = h3 / r1  #  velocity at periapsis for the Hohmann

print 'The velocity at periapsis for the Hohmann transfer is ' + str(vp3)

va3 = h3 / r2

print 'The velocity at apoapsis for the Hohmann transfer is ' + str(va3)

deltavA = np.absolute(vp3 - v1)

print 'The Delta v requirement at A is ' + str(deltavA)

T2 = 2.0 * np.pi / np.sqrt(mue) * np.sqrt(r2 ** 3)

phaseS2 = tofsec / T2 * 360 - 180.0

print 'S2 lags S1 by ' + str(phaseS2) + ' degrees'

Tphase = T2 * (1.0 - (phaseS2 / 360.0) / 3.0)
#  S1 one needs an orbit that is longer and will execute phaseS2 in
#  three revolutions

print 'The period of the phase orbit is then ' + str(Tphase)

aphase = (Tphase * np.sqrt(mue) / (2.0 * np.pi)) ** (2.0 / 3.0)
#  The semi-major axis for the phase

print 'The semi-major axis for the phasing orbit is ' + str(aphase)

raphase = 2.0 * aphase - r2

print 'The apoapsis of the new orbit is ' + str(raphase)

ecc = (raphase - r2) / (raphase + r2)

print 'The eccentricity of the new orbit is ' + str(ecc)

hphase = np.sqrt(2.0 * mue) * np.sqrt((raphase * r2) / (raphase + r2))

print 'The angular momentum of the phasing orbit at periapsis is ' + str(hphase)

vphaseC = hphase / r2  #  velocity at periapsis

print 'The velocity of phasing orbit at C is ' + str(vphaseC)

deltavC = np.absolute(vphaseC - va3)

print 'The Delta v at C is ' + str(deltavC)

deltav2 = np.absolute(v2 - vphaseC)

print 'The final delta V required to slow down at C is ' + str(deltav2)

DeltaV = deltav2 + deltavC + deltavA

print 'The total Delta V is ' + str(DeltaV)
