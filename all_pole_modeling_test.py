# -*- coding: utf-8 -*-
import math
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import numpy.linalg
import scipy.signal
import os
from pylab import *

pi = math.pi

# the unit circle (for plotting)
x_uc = np.linspace(-1,1,100)
ypos_uc = (1-x_uc**2)**0.5
yneg_uc = -ypos_uc

# target signal
n = np.arange(100)
#x = 0.5**n
x = np.sin(2*pi/10*n)
x = np.concatenate((x[3:],x[:3])) # shift so that x[0] is non-zero

# deterministic autocorrelation
rxx = []
for w in range(2*len(x)-1):
    val = 0
    if w<5:
        print 'w = ',w
    for d in range(len(x)):
        if -1<99-w+d<100:
            if w<5:
                print '***  ',d,99-w+d
            val += x[d]*x[99-w+d]
    rxx.append(val)

# model order
#p = 1 # Ediff =  51.5207233409
#W = rxx[99]
#D = rxx[100]
p = 2 # Ediff =  57.5531258058
W = np.array([[rxx[99],rxx[100]],[rxx[100],rxx[99]]])
D = np.array([[rxx[100]],[rxx[101]]])
#p = 3 # Ediff = 63.76142098
#W = np.array([[rxx[99],rxx[100],rxx[101]],[rxx[100],rxx[99],rxx[109]],[rxx[191],rxx[109],rxx[99]]])
#D = np.array([[rxx[100]],[rxx[101]],[rxx[102]]])
#p = 4 # Ediff =  107.494904876
#W = np.array([[rxx[99],rxx[100],rxx[101],rxx[102]],[rxx[100],rxx[99],rxx[100],rxx[101]],[rxx[101],rxx[100],rxx[99],rxx[100]],
#             [rxx[102],rxx[101],rxx[100],rxx[99]]])
#D = np.array([[rxx[100]],[rxx[101]],[rxx[102]],[rxx[103]]])

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
print 'norms of zeros = ',np.abs(zeros)
print 'angles of zeros = ',np.angle(zeros)
print '*** compare with 2*pi/10 = ',2*pi/10
real = [elt.real for elt in zeros]
imag = [elt.imag for elt in zeros]

# form impulse response of all-pole filter
# (a.k.a. all-pole model of target signal)
num = 1.
den = a # [1.,0.]
dt = 1.
t,h = scipy.signal.dimpulse((num,den,dt),n=len(x)+1)
h = list(h[0])
h = [float(elt) for elt in h][1:] # for some reason h[0] is always 0
print 'h[:5] = ',h[:5]
print 'h.index(max(h)) = ',h.index(max(h))
# compare model to target signal
diff = x-h
Ediff = sum(diff**2)
print 'Ediff = ',Ediff

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

fig3 = plt.figure()
ax4 = fig3.add_subplot(211)
ax4.stem(h,label='h[n]')
ax4.set_xlim(left=-1,right=len(h))
ax4.set_ylim(top=max(h)*1.2)
ax4.legend(loc='best',frameon=False)
ax5 = fig3.add_subplot(212)
ax5.stem(x,label='x[n]')
ax5.set_xlim(left=-1,right=len(x))
ax5.set_ylim(top=max(x)*1.2)
ax5.legend(loc='best',frameon=False)

plt.show()
