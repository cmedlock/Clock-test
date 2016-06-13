# -*- coding: utf-8 -*-
import math
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import numpy.linalg
import numpy.fft
import os
from pylab import *

import clock_test as ct
ct = reload(ct)

path = '/Users/cmedlock/Documents/DSP_UROP/simulated data/'
dirs = os.listdir(path)

# path to DFS coefficient figures
if not os.path.exists(path+'compare_circles_ellipses'):
    os.makedirs(path+'compare_circles_ellipses')

pi = math.pi

# fraction of energy contained in fundamental frequency
Epeak_x,Epeak_y = [],[]

axis_ratios = np.linspace(1,1.2,100)
for ratio in axis_ratios:

    x_eqdist = np.loadtxt(path+'norm velocity data/perfect_ellipse_x_eqdist_ratio_'+str(ratio)[:5]+'.txt')
    y_eqdist = np.loadtxt(path+'norm velocity data/perfect_ellipse_y_eqdist_ratio_'+str(ratio)[:5]+'.txt')

    # DFS coefficients
    dft_size = len(x_eqdist)
    dftx,dfty = np.fft.fft(x_eqdist,n=dft_size),np.fft.fft(y_eqdist,n=dft_size)
    k = np.arange(dft_size)

    # k_near_pi is the k value for which w_k = 2*pi*k/N is closest to,
    # but not larger than, pi
    k_near_pi = 0
    if dft_size%2==0:
        k_near_pi = dft_size/2+1
    else:
        k_near_pi = math.ceil(dft_size/2)

    # percent energy in peak
    Ex,Ey = np.abs(dftx)**2,np.abs(dfty)**2
    Ex_total,Ey_total = sum(Ex),sum(Ey)
    Ex_peak,Ey_peak = 2*Ex[1]/Ex_total,2*Ey[1]/Ey_total
    Epeak_x.append(Ex_peak)
    Epeak_y.append(Ey_peak)

# compare fraction of energy in fundamental frequency
plt.close('all')
fig = plt.figure()
fig.subplots_adjust(hspace=0.6,left=0.15,bottom=0.15)
ax = fig.add_subplot(111)
ax.clear()

Epeak_x_low_eccentricity,Epeak_y_low_eccentricity = [],[]
Epeak_x_high_eccentricity,Epeak_y_high_eccentricity = [],[]
for w in range(len(axis_ratios)):
    if axis_ratios[w]<1.05:
        Epeak_x_low_eccentricity.append(Epeak_x[w])
        Epeak_y_low_eccentricity.append(Epeak_y[w])
    else:
        Epeak_x_high_eccentricity.append(Epeak_x[w])
        Epeak_y_high_eccentricity.append(Epeak_y[w])

ax.scatter(Epeak_x_low_eccentricity,Epeak_y_low_eccentricity,color='green',marker='o',alpha=0.5,label='(major axis)/(minor axis) < 1.05')
ax.scatter(Epeak_x_high_eccentricity,Epeak_y_high_eccentricity,color='black',marker='o',alpha=0.5,label='(major axis)/(minor axis) > 1.05')
#ax.scatter(axis_ratios[:len(Epeak_x_low_eccentricity)],Epeak_y_low_eccentricity,color='green',marker='o',alpha=0.5,label='(major axis)/(minor axis) < 1.05')
#ax.scatter(axis_ratios[len(Epeak_x_low_eccentricity):],Epeak_y_high_eccentricity,color='black',marker='o',alpha=0.5,label='(major axis)/(minor axis) > 1.05')
#ax.set_xlim(min(Epeak_x),max(Epeak_x))
#ax.set_ylim(min(Epeak_y),max(Epeak_y))
ax.set_xlim(0.9994,1.000)
ax.set_ylim(0.9992,1.000)
y_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
ax.yaxis.set_major_formatter(y_formatter)
ax.xaxis.set_major_formatter(y_formatter)
ax.set_xlabel(r'$E_{peak}^x/E_{total}$',fontsize=20)
ax.set_ylabel(r'$E_{peak}^y/E_{total}$',fontsize=20)
ax.legend(loc='best')
fig.savefig(path+'compare_circles_ellipses/compare_Epeak_x_and_y.png')
