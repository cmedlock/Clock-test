# switch from rectangular to polar coordinates

import math
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import numpy.fft
import os
from pylab import *

# get coordinates and timestamps
n = np.arange(250)
x = np.cos(2*np.pi/250*n)
y = np.sin(2*np.pi/250*n)

# find COM
x_com = np.mean(x)
y_com = np.mean(y)

# get r and theta
r,theta = [],[]
for w in range(len(x)):
    dx,dy = x[w]-x_com,y[w]-y_com
    dist = sqrt(dx**2+dy**2)
    angle = math.atan2(dy,dx)
    r.append(dist)
    theta.append(angle)

# plot
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(theta,r)
ax.set_xlabel('theta')
ax.set_ylabel('r')
plt.show()
