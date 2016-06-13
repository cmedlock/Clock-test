# -*- coding: utf-8 -*-
import math
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import numpy.fft
import os
from pylab import *

import clock_test as ct
ct = reload(ct)

path = '/Users/cmedlock/Documents/DSP_UROP/simulated data/'
dirs = os.listdir(path)

# path to figures for each individual file
if not os.path.exists(path+'compare circles ellipses'):
    os.makedirs(path+'compare circles ellipses')

pi = math.pi

# variance of short time frequency distribution
freq_var_x,freq_var_y = [],[]

axis_ratios = np.linspace(1,1.2,100)
for ratio in axis_ratios:

    x_eqdist = np.loadtxt(path+'norm velocity data/perfect_ellipse_x_eqdist_ratio_'+str(ratio)[:5]+'.txt')
    y_eqdist = np.loadtxt(path+'norm velocity data/perfect_ellipse_y_eqdist_ratio_'+str(ratio)[:5]+'.txt')

    # downsample by a factor of 3
    x_downsampled,y_downsampled = [],[]
    for w in range(len(x_eqdist)):
        if w%3==0:
            x_downsampled.append(x_eqdist[w])
            y_downsampled.append(y_eqdist[w])
    #x_eqdist,y_eqdist = np.array(x_downsampled),np.array(y_downsampled)

    # extract frequency from each set of 3 points using the formula
    # w[n] = 2*cos(omega)*w[n-1] - w[n-2]
    x_freq,y_freq = [],[]
    for w in range(len(x_eqdist)):
        frequency_x = (x_eqdist[w]+x_eqdist[w-2])/(2*x_eqdist[w-1])
        frequency_y = (y_eqdist[w]+y_eqdist[w-2])/(2*y_eqdist[w-1])
        x_freq.append(frequency_x)
        y_freq.append(frequency_y)
    freq_var_x.append(var(x_freq))
    freq_var_y.append(var(y_freq))
    
# compare short time frequency variance
plt.close('all')
fig = plt.figure()
fig.subplots_adjust(hspace=0.6,left=0.15,bottom=0.15)
ax = fig.add_subplot(111)
ax.clear()

freq_var_x_low_eccentricity,freq_var_y_low_eccentricity = [],[]
freq_var_x_high_eccentricity,freq_var_y_high_eccentricity = [],[]
for w in range(len(axis_ratios)):
    if axis_ratios[w]<1.05:
        freq_var_x_low_eccentricity.append(freq_var_x[w])
        freq_var_y_low_eccentricity.append(freq_var_y[w])
    else:
        freq_var_x_high_eccentricity.append(freq_var_x[w])
        freq_var_y_high_eccentricity.append(freq_var_y[w])
"""
ax.scatter(freq_var_x_low_eccentricity,freq_var_y_low_eccentricity,color='green',marker='o',alpha=0.5,label='(major axis)/(minor axis) < 1.05')
ax.scatter(freq_var_x_high_eccentricity,freq_var_y_high_eccentricity,color='black',marker='o',alpha=0.5,label='(major axis)/(minor axis) > 1.05')
#ax.scatter(axis_ratios[:len(freq_var_x_low_eccentricity)],freq_var_y_low_eccentricity,color='green',marker='o',alpha=0.5,label='(major axis)/(minor axis) < 1.05')
#ax.scatter(axis_ratios[len(freq_var_x_low_eccentricity):],freq_var_y_high_eccentricity,color='black',marker='o',alpha=0.5,label='(major axis)/(minor axis) > 1.05')
ax.set_xlim(2*10**-8,10**-7)
#ax.set_ylim(0.9*min(freq_var_y),1.1*max(freq_var_y))
ax.set_xlabel(r'var$(f(n))$ for $x[n]$',fontsize=20)
ax.set_ylabel(r'var$(f(n))$ for $y[n]$',fontsize=20)
j=1
for v1 in ax.get_xticklabels():
    v1.set_fontsize(17)
    if j==1:
        v1.set_visible(False)
        j=0
for v2 in ax.get_yticklabels():
    v2.set_fontsize(17)
ax.legend(loc='best')
fig.savefig(path+'compare circles ellipses/compare_freq_var_x_and_y.png')
"""
ax.clear()
ax.scatter(axis_ratios[:len(freq_var_x_low_eccentricity)],freq_var_x_low_eccentricity,color='green',marker='o',alpha=0.5,label='(major axis)/(minor axis) < 1.05')
ax.scatter(axis_ratios[len(freq_var_x_low_eccentricity):],freq_var_x_high_eccentricity,color='black',marker='o',alpha=0.5,label='(major axis)/(minor axis) > 1.05')
#ax.set_xlim(2*10**-8,10**-7)
ax.set_ylim(2*10**-8,10**-7)
ax.set_xlabel('(major axis)/(minor axis)',fontsize=20)
ax.set_ylabel(r'var$(f(n))$ for $x[n]$',fontsize=20)
j=1
for v1 in ax.get_xticklabels():
    v1.set_fontsize(17)
    if j==1:
        v1.set_visible(False)
        j=0
for v2 in ax.get_yticklabels():
    v2.set_fontsize(17)
ax.legend(loc='best')
fig.savefig(path+'compare circles ellipses/compare_freq_var_x.png')

ax.clear()
ax.scatter(axis_ratios[:len(freq_var_y_low_eccentricity)],freq_var_y_low_eccentricity,color='green',marker='o',alpha=0.5,label='(major axis)/(minor axis) < 1.05')
ax.scatter(axis_ratios[len(freq_var_y_low_eccentricity):],freq_var_y_high_eccentricity,color='black',marker='o',alpha=0.5,label='(major axis)/(minor axis) > 1.05')
#ax.set_xlim(2*10**-8,10**-7)
ax.set_ylim(10,25)
ax.set_xlabel('(major axis)/(minor axis)',fontsize=20)
ax.set_ylabel(r'var$(f(n))$ for $y[n]$',fontsize=20)
j=1
for v1 in ax.get_xticklabels():
    v1.set_fontsize(17)
    if j==1:
        v1.set_visible(False)
        j=0
for v2 in ax.get_yticklabels():
    v2.set_fontsize(17)
ax.legend(loc='best')
fig.savefig(path+'compare circles ellipses/compare_freq_var_y.png')
