#!/usr/bin/env ipython
#  This program solves example 17 of Orbital 230 Notes.

import numpy as np

me = 5.974e24  #  mass of the Earth
ms = 1.989e30  #  mass of the Sun
Re = 1.496e8  #  radius of the Earth's orbit

rsoi = Re * (me / ms) ** (2.0 / 5.0)  #  SOI Earth

print 'The sphere of influence of the Earth is ' + str(rsoi)
