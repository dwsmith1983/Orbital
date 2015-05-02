#!/usr/bin/env ipython
#  This program solves example 11 of Orbital Notes 230.
#  Circular Hohmann transfer Earth to Mars

import numpy as np

aE = 149598261.0  #  semi-major axis Earth
aM = 227939100.0  #  semi-major axis Mars
musun = 132712440018  #  grav param of the Sun

a = .5 * (aE + aM)  #  semi-major axis for the transfer ellipse

print 'The semi-major axis for the transfer ellipse is ' + str(a)

tofs = np.pi / np.sqrt(musun) * np.sqrt(a ** 3)
#  time of flight in secs

tofd = tofs / (3600.0 * 24.0)
#  time of flight in days

print 'The time of flight in seconds is ' + str(tofs)
print 'The time of flight in days is ' + str(tofd)

TMars = 2 * np.pi * np.sqrt(aM ** 3) / np.sqrt(musun)
#  period of Mars

alpha = 180 - (tofs / TMars * 360)
#  The lead angle for the Earth to Mars transfer

print 'The lead angle of Mars is ' + str(alpha)
