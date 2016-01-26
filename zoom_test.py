# -*- coding: utf-8 -*-
# sanity check for 'zoom' transform

import math
import cmath
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
        
path = '/Users/cmedlock/Documents/DSP_UROP/all_data/'
dirs = os.listdir(path)

if not os.path.exists(path+'figs_raw'):
    os.makedirs(path+'figs_raw')

n = np.arange(6)
x = np.array([1]*6)

pi = math.pi
dftx = np.fft.fft(x,n=500)
k = np.arange(500)
# check
dftx_check = [sum(x)]
for w in range(1,6):
    X_of_k = cmath.exp(-1j*5*pi*w/6)*math.sin(pi*w)/math.sin(pi*w/6)
    if abs(X_of_k)<10**-12:
        X_of_k = 0
    dftx_check.append(X_of_k)
dftx_check = np.array(dftx_check)
print dftx==dftx_check
# another check
X_ejw = np.exp(-1j*pi*k/12)*np.sin(pi*k/10)/np.sin(pi*k/60)

# zoom in on the frequencies between ±(pi/3)
# convolve (truncated) sinc with x[n]
omega_c = pi/3
length_of_sinc = 201
x1 = []
for n in range(len(x)+length_of_sinc-1):
    x1_of_n = 0
    for m in range(len(x)):
        x1_of_n += x[m]*sinc(omega_c,n-100-m,length_of_sinc)
    x1.append(x1_of_n)
print 'check: ',x[0]*sinc(omega_c,100,length_of_sinc),' = ',x1[0],'?'
x1_of_0 = 0
for d in range(6):
    x1_of_0 += sinc(omega_c,d,length_of_sinc)
print 'check: ',x1_of_0,' = ',x1[100],'?'
# check x1[n]
dftx1 = np.fft.fft(x1,n=500)
k1 = np.arange(500)
print 'check: ',sum(x1),' = ',3*sum(x),'?'

# check filter
h1 = []
for w in range(-100,101):
    h1.append(sinc(omega_c,w,length_of_sinc))
dfth1 = np.fft.fft(h1,n=500)
kh1 = np.arange(500)

# downsample by 3
for w in range(len(x1)-1):
    x1 = np.insert(x1,w*3+1,[0]*2)

# subtract off mean
#x1_zm = x1-mean(x1)

# redo the N-point DFT
dftx_zoom = np.fft.fft(x1,n=len(x))
for w in range(len(dftx_zoom)):
    if np.abs(dftx_zoom[w])<10**-10:
        dftx_zoom[w] = 0
# check
for d in range(len(dftx_zoom)):
    if d%3==0:
        print 'd = ',d,': ',np.abs(dftx[d/3]),np.abs(dftx_zoom[d])

plt.close('all')
fig = plt.figure()
ax1 = fig.add_subplot(311)
ax1.stem(k,np.abs(dftx))
ax1.set_xlim(left=-1,right=max(k)+1)
ax1.set_ylim(bottom=-1,top=max(abs(dftx))+1)
ax2 = fig.add_subplot(312)
ax2.stem(kh1,np.abs(dfth1))
ax2.set_xlim(left=-1,right=max(kh1)+1)
ax2.set_ylim(bottom=-1,top=max(np.abs(dfth1)+1))
ax3 = fig.add_subplot(313)
ax3.stem(k1,np.abs(dftx1))
ax3.set_xlim(left=-1,right=max(k1)+1)
ax3.set_ylim(bottom=-1,top=max(np.abs(dftx1)+1))
plt.show()
