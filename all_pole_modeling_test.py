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
n = np.arange(100)
x = 0.5**n

# deterministic autocorrelation
rxx = []
for w in range(2*len(x)-1):
    val = 0
    for d in range(len(x)):
        if -1<99-w+d<100:
            val += x[d]*x[99-w+d]
    rxx.append(val)

# model order
#p = 1 # fraction of E_g in g[0] =  1.0
#W = rxx[99]
#D = rxx[100]
#p = 2 # fraction of E_g in g[0] =  1.0
#W = np.array([[rxx[99],rxx[100]],[rxx[100],rxx[99]]])
#D = np.array([[rxx[100]],[rxx[101]]])
#p = 3 # fraction of E_g in g[0] =  0.977961572123
#W = np.array([[rxx[99],rxx[100],rxx[101]],[rxx[100],rxx[99],rxx[109]],[rxx[191],rxx[109],rxx[99]]])
#D = np.array([[rxx[100]],[rxx[101]],[rxx[102]]])
p = 4 # fraction of E_g in g[0] =  1.0
W = np.array([[rxx[99],rxx[100],rxx[101],rxx[102]],[rxx[100],rxx[99],rxx[100],rxx[101]],[rxx[101],rxx[100],rxx[99],rxx[100]],
             [rxx[102],rxx[101],rxx[100],rxx[99]]])
D = np.array([[rxx[100]],[rxx[101]],[rxx[102]],[rxx[103]]])

if p>1:
    W_inv = np.linalg.inv(W)
    ak = np.dot(W_inv,D) # linear prediction coefficients
else:
    W_inv = 1/W
    ak = [W_inv*D]

# form impulse response of linear prediction filter
a = [1]
for w in range(p):
    a.append(-float(ak[w]))

# convolve with x[n]
g = []
for w in range(len(x)+p):
    val = 0
    for d in range(len(x)):
        if -1<w-d<p+1:
            val += x[d]*a[w-d]
    g.append(val)
g = np.array(g)
print 'fraction of E_g in g[0] = ',g[0]**2/sum(g**2)

# plot zeros of linear prediction filter
zeros = 1/np.roots(a[::-1])
real = [elt.real for elt in zeros]
imag = [elt.imag for elt in zeros]

plt.close('all')

fig1 = plt.figure()
ax1 = fig1.add_subplot(211)
ax1.stem(g,label='g[n]')
ax1.set_xlim(left=-1,right=len(g))
ax1.set_ylim(top=max(g)*1.2)
ax1.legend(loc='best',frameon=False)
ax2 = fig1.add_subplot(212)
ax2.stem(x,label='x[n]')
ax2.set_xlim(left=-1,right=len(x))
ax2.set_ylim(top=max(x)*1.2)
ax2.legend(loc='best',frameon=False)

fig2 = plt.figure()
ax3 = fig2.add_subplot(111)
ax3.scatter(real,imag)
ax3.plot(x_uc,ypos_uc,'darkgreen')
ax3.plot(x_uc,yneg_uc,'darkgreen')
ax3.axvline(x=0, color='b')
ax3.axhline(y=0, color='b')
ax3.set_xlabel('real')
ax3.set_ylabel('imag')
plt.axis('equal')

plt.show()
