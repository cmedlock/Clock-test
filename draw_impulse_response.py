# -*- coding: utf-8 -*-
# 1.normalize for non-constant velocity
# 2.estimate ideal underlying sinusoids
# 3.compute correlation and mean squared difference
# between ideal sinusoids and actual drawing
import math
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import numpy.linalg
import scipy.signal
import os
from pylab import *

import clock_test as ct
ct = reload(ct)

pi = math.pi

# the unit circle (for plotting)
x_uc = np.linspace(-1,1,100)
ypos_uc = (1-x_uc**2)**0.5
yneg_uc = -ypos_uc

# a test filter
num = [-1./4.,0.,1.]
den = 1.
tout,yout = scipy.signal.dimpulse((num,den,[1.]))
print 'tout = ',tout
print 'yout = ',yout

plt.close('all')

fig1 = plt.figure()
ax1 = fig1.add_subplot(211)
ax1.stem(g,label='g[n]')
ax1.set_xlim(left=-1,right=12)
ax1.legend(loc='best',frameon=False)
ax2 = fig1.add_subplot(212)
ax2.stem(x,label='x[n]')
ax2.set_xlim(left=-1,right=12)
ax2.legend(loc='best',frameon=False)

fig2 = plt.figure()
ax3 = fig2.add_subplot(111)
ax3.scatter(real,imag)
ax3.plot(x_uc,ypos_uc,'darkgreen')
ax3.plot(x_uc,yneg_uc,'darkgreen')
ax3.set_xlabel('real')
ax3.set_ylabel('imag')
plt.axis('equal')
plt.show()