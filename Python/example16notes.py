#!\usr\bin\env ipython
#  This program solve problem 16 of Orbital Notes 230

import numpy as np

musun = 1.32712e11  #  grav param
aE = 1.496e8  #  semi-major axis of Earth
aV = 1.082e8  #  semi-major axis of Venus
phii = 190.88275 * np.pi / 180.0  #  initial lead angle in radians
th = np.pi / np.sqrt(musun) * np.sqrt(((aE + aV) / 2) ** 3)
thdays = th / (3600.0 * 24.0)

print 'The time of the Hohmann transfer is ' + str(th) + ' sec'

print 'The time in days is ' + str(thdays)

TV = 224.7 * 24.0 * 3600.0  #  the period of Venus in secs
TE = 365.25 * 24.0 * 3600.0  #  the period of Earth in secs
nV = 2.0 * np.pi / TV  #  mean motion of Venus
nE = 2.0 * np.pi / TE  #  mean motion of Earth
phi0 = 2.0 * np.pi + np.pi - nV * th
#  the phase angle at launch from Venus to Earth
phi0degrees = phi0 / np.pi * 180.0

print 'The phase angle launch is ' + str(phi0degrees) + ' degrees'

tlaunch = (phi0 - phii) / (nV - nE)  #  the time to launch in sec
tlaunchdays = tlaunch / (3600.0 * 24.0)

print 'The time to launch is ' + str(tlaunch) + ' sec'

print 'The launch time in days is ' + str(tlaunchdays)

nufE = nE * th  #  true anomaly in radians

print 'The true anomaly of Earth is ' + str(nufE / np.pi * 180.0) + ' degrees'

phif = np.pi - nufE

print 'The final phase angle is ' + str(phif / np.pi * 180.0) + ' degrees'

twait = (2 * np.pi - 2 * phif) / (nV - nE)  #  wait time in sec
twaitdays = twait / (24.0 * 3600.0)

print 'The wait time is ' + str(twaitdays) + ' days'

totaltime = tlaunch + 2 * th + twait  #  total time sec
totaltimedays = tlaunchdays + 2 * thdays + twaitdays
totaltimeyears = totaltime / (24.0 * 3600.0 * 365.25)

print 'The total launch time is ' + str(totaltime) + ' sec'

print 'The total time in days is ' + str(totaltimedays)

print 'The total time in years is ' + str(totaltimeyears)
