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
bigN = 69
print 'N = ',bigN

pi = math.pi
dftx = np.fft.fft(x,n=bigN)
k = np.arange(bigN)
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
# check some individual values
print 'check: ',x[0]*sinc(omega_c,100,length_of_sinc),' = ',x1[0],'?'
x1_of_0 = 0
for d in range(6):
    x1_of_0 += sinc(omega_c,d,length_of_sinc)
print 'check: ',x1_of_0,' = ',x1[100],'?'

# check dft of x1[n]
dftx1 = np.fft.fft(x1,n=bigN)
k1 = np.arange(bigN)
print 'check: ',sum(x1),' = ',3*sum(x),'?'

# check filter
h1 = []
for w in range(-100,101):
    h1.append(sinc(omega_c,w,length_of_sinc))
dfth1 = np.fft.fft(h1,n=bigN)
kh1 = np.arange(bigN)

# downsample by 3 (d = decimated)
x1_d = []
for w in range(len(x1)):
    if w%3==0:
        x1_d.append(x1[w])

# redo the N-point DFT
dftx_zoom = np.fft.fft(x1_d,bigN)
print 'check: ',sum(x1_d),' = ',np.abs(dftx_zoom[0]),'?'
print 'or is it: ',sum(x1_d[:bigN]),' = ',np.abs(dftx_zoom[0]),'?'
for w in range(len(dftx_zoom)):
    if np.abs(dftx_zoom[w])<10**-10:
        dftx_zoom[w] = 0
# check
for d in range(len(dftx_zoom)):
    if d%3==0 and d < 30:
        print 'd = ',d,': ',np.abs(dftx[d/3]),np.abs(dftx_zoom[d]),' -> ',abs(np.abs(dftx[d/3])-np.abs(dftx_zoom[d]))
# for bigN = 25
dftx_zoom_at0_check = 0

plt.close('all')
fig = plt.figure()
ax1 = fig.add_subplot(411)
ax1.stem(k,np.abs(dftx))
ax1.set_xlim(left=-1,right=max(k)+1)
ax1.set_ylim(bottom=-1,top=max(abs(dftx))+1)
ax2 = fig.add_subplot(412)
ax2.stem(kh1,np.abs(dfth1))
ax2.set_xlim(left=-1,right=max(kh1)+1)
ax2.set_ylim(bottom=-1,top=max(np.abs(dfth1))+1)
ax3 = fig.add_subplot(413)
ax3.stem(x1_d[:bigN])
ax3.set_xlim(left=-1,right=len(x1_d[:bigN])+1)
ax3.set_ylim(bottom=-1,top=max(x1_d[:bigN])+1)
ax4 = fig.add_subplot(414)
ax4.stem(k,np.abs(dftx_zoom))
ax4.set_xlim(left=-1,right=max(k)+1)
ax4.set_ylim(bottom=-1,top=max(np.abs(dftx_zoom))+1)
plt.show()
