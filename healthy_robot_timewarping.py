# for sanity checks on normalizing for non-constant velocity

import math
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import numpy.fft
import os
from pylab import *

import clock_test as ct
ct = reload(ct)

def sinc(omega_c,n,length_of_sinc):
    if n==0:
        return 1
    elif abs(n)>(length_of_sinc-1)/2:
        return 0
    else:
        return math.sin(omega_c*n)/(omega_c*n)

# get coordinates and timestamps
n = np.arange(141)
x = np.cos(1./2.*(2.*np.pi*n/250.)**2)
y = np.sin(1./2.*(2.*np.pi*n/250.)**2)

# velocity/distance (T = 1)
dists = []
for w in range(1,len(x)):
    dx,dy = x[w]-x[w-1],y[w]-y[w-1]
    dist = math.sqrt(dx**2+dy**2)
    dists.append(dist)
dist_avg = mean(dists)
print 'average distance between points is ',dist_avg

# want to get 141 evenly-spaced points along the curve

# generate a much longer array with 3 linearly-interpolated 
# points between the actual data points
x_interp,y_interp = [],[]
for w in range(len(x)-1):
    x_interp.append(x[w])
    y_interp.append(y[w])
    dx,dy = x[w+1]-x[w],y[w+1]-y[w]
    dist = math.sqrt(dx**2+dy**2)
    for r in range(1,10):
        x_new = x[w]+r*dx/10
        y_new = y[w]+r*dy/10
        x_interp.append(x_new)
        y_interp.append(y_new)
x_interp.append(x[-1])
y_interp.append(y[-1])
print 'average distance between interpolated points is ',dist_avg/10

# start from the first point and find the ones that are 
# approximately a distance dist_avg from each other
x_eqdist,y_eqdist = [x_interp[0]],[y_interp[0]]
idx = 0 
for k in range(len(x)):
    dist_total = 0
    for j in range(idx,len(x_interp)-1):
        dx,dy = x_interp[j+1]-x_interp[j],y_interp[j+1]-y_interp[j]
        dist_total += math.sqrt(dx**2+dy**2)
        if abs(dist_total-dist_avg)<0.01:
            idx = j+1
            break
    x_eqdist.append(x_interp[idx])
    y_eqdist.append(y_interp[idx])
print len(x_eqdist),len(y_eqdist)
plt.close('all')

fig_xy = plt.figure()
ax_xy = fig_xy.add_subplot(111)
ax_xy.plot(x,y)
ax_xy.set_xlabel('x')
ax_xy.set_ylabel('y')

fig = plt.figure()
ax1 = fig.add_subplot(311)
ax1.plot(x,label='x')
ax1.plot(y,label='y')
ax1.legend(loc='best',frameon=False)
ax1.set_ylabel('x,y',fontsize=20)
ax1.set_ylim(bottom=min(x)-0.2,top=max(x)+0.2)
ax2 = fig.add_subplot(312)
ax2.plot(dists)
ax2.set_ylabel('v',fontsize=20)
ax3 = fig.add_subplot(313)
ax3.plot(x_eqdist,label='x_eqdist')
ax3.plot(y_eqdist,label='y_eqdist')
ax3.legend(loc='best',frameon=False)
ax3.set_ylabel('x_eqdist,\ny_eqdist',fontsize=15)
ax3.set_xlim(right=len(x_eqdist))
plt.show()
