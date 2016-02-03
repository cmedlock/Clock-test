# -*- coding: utf-8 -*-
import math
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import numpy.linalg
import os
from pylab import *

pi = math.pi

# the unit circle (for plotting)
x_uc = np.linspace(-1,1,100)
ypos_uc = (1-x_uc**2)**0.5
yneg_uc = -ypos_uc

# target signal
n = np.arange(10)
x = np.sin(2*pi/10*n)
x = np.concatenate((x[3:],x[:3])) # shift so that x[0] is non-zero

# deterministic autocorrelation
rxx = []
for w in range(19):
    val = 0
    for d in range(10):
        if -1<9-w+d<10:
            val += x[d]*x[9-w+d]
    rxx.append(val)

# model order
#p = 2
#W = np.array([[rxx[9],rxx[10]],[rxx[10],rxx[9]]])
#D = np.array([[rxx[10]],[rxx[11]]])
#p = 3
#W = np.array([[rxx[9],rxx[10],rxx[11]],[rxx[10],rxx[9],rxx[10]],[rxx[11],rxx[10],rxx[9]]])
#D = np.array([[rxx[10]],[rxx[11]],[rxx[12]]])
p = 4
W = np.array([[rxx[9],rxx[10],rxx[11],rxx[12]],[rxx[10],rxx[9],rxx[10],rxx[11]],[rxx[11],rxx[10],rxx[9],rxx[10]],[rxx[12],rxx[11],rxx[10],rxx[9]]])
D = np.array([[rxx[10]],[rxx[11]],[rxx[12]],[rxx[13]]])

W_inv = np.linalg.inv(W)
ak = np.dot(W_inv,D) # linear prediction coefficients

# form impulse response of linear prediction filter
a = [1]
for w in range(p):
    a.append(-float(ak[w]))

# convolve with x[n]
g = []
for w in range(10): # only go up to 10, because that's as far as we're interested in predicting/modeling
    val = 0
    for d in range(10):
        if -1<w-d<3:
            val += x[d]*a[w-d]
    g.append(val)
g = np.array(g)
print 'fraction of E_g in g[0] = ',g[0]**2/sum(g**2)

# plot zeros of linear prediction filter
zeros = np.roots(a[::-1])
real = [elt.real for elt in zeros]
imag = [elt.imag for elt in zeros]

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
