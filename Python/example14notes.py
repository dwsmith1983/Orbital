#!/usr/bin/env ipython
#  This program solves example 14 from Orbital Notes 230

import numpy as np

rLEO = 6371.0 + 278.0
rGEO = 6371.0 + 35786.0
mue = 398600.0
angle = 28.0 / 180.0 * np.pi  #  in radians

vLEO = np.sqrt(mue / rLEO)
vGEO = np.sqrt(mue / rGEO)

a = .5 * (rLEO + rGEO)  #  semi-major axis

tofsec = np.pi / np.sqrt(mue) * np.sqrt(a ** 3)
tofhours = tofsec / 3600.0

print 'The total time of flight is ' + str(tofhours) + ' hours'

h = np.sqrt(2.0 * mue) * np.sqrt((rLEO * rGEO) / (rLEO + rGEO))

vp = h / rLEO

deltav1 = np.absolute(vp - vLEO)

va = h / rGEO

deltav2 = np.sqrt(vGEO ** 2 + va ** 2 - 2 * vGEO * va * np.cos(angle))

print 'The total Delta v is ' + str(deltav1 + deltav2)
