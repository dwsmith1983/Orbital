#!/usr/bin/env ipython
#  This program solves the bi-elliptical Earth transfer for
#  230 Notes.

import numpy as np

re = 6371.0  #  radius of Earth in KM
r1 = 300.0 + re
r4 = 3000.0 + re
ecc2 = 0.3
mue = 398600.0  #  grav param of Earth

v1 = np.sqrt(mue / r1)  #  speed of satellite

print 'The speed of the satellite at 300 km is ' + str(v1)

#  angular moment of 2
h2 = np.sqrt(mue * r1 * (1 + ecc2))

print 'The angular momentum of the 1st transfer ellipse is ' + str(h2)

v2p = h2 / r1  #  speed transfer ellipse 2

print 'The speed needed to execute the 1st transfer ellipse is ' + str(v2p)

deltavA = np.absolute(v2p - v1)

print 'Delta vA is ' + str(deltavA)

ra2 = h2 ** 2 / (mue * (1 - ecc2))  #  radius at apoapsis of orbit 2

print 'The radius at apoapsis of orbit 2 is ' + str(ra2)

v2a = h2 / ra2  #  speed at apoapsis

print 'The speed of the satellite at apoapsis is ' + str(v2a)

ecc3 = (ra2 - r4) / (ra2 + r4)  #  eccentricity of the 2nd ellipse

print 'The eccentricity of the 2nd transfer ellipse is ' + str(ecc3)

h3 = np.sqrt(mue * r4 * (1 + ecc3))

print 'The angular momentum at B for the 2nd ellipse is ' + str(h3)

v3a = h3 / ra2

print 'The speed at apoapsis of ellipse 2 is ' + str(v3a)

deltavB = np.absolute(v3a - v2a)

print 'Delta vB is ' + str(deltavB)

v3p = h3 / r4

print 'The speed at arrival at C is ' + str(v3p)

v4 = np.sqrt(mue / r4)

print 'The speed of the satellite at 3000 km is ' + str(v4)

deltavC = np.absolute(v4 - v3p)

print 'Delta vC is ' + str(deltavC)

deltav = deltavA + deltavB + deltavC

print 'The total delta v is ' + str(deltav)

a2 = ra2 / (1 + ecc2)  #  semi-major axis 2

a3 = ra2 / (1 + ecc3)  #  semi-major axis 3

tofsec = np.pi / np.sqrt(mue) * (np.sqrt(a2 ** 3) + np.sqrt(a3 ** 3))

print 'The time of flight is ' + str(tofsec) + ' sec'

tofhours = tofsec / 3600.0

print 'The time of flight is ' + str(tofhours) + ' hours'
