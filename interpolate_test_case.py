# -*- coding: utf-8 -*-
import math
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import numpy.fft
import os
from pylab import *

path = '/Users/cmedlock/Documents/DSP_UROP/all_data/'

if not os.path.exists(path+'figs_raw'):
    os.makedirs(path+'figs_raw')

# length of truncated sinc (longer = more accurate interpolation)
length = 20
missing_points = [5,15]

# original signal, appropriately bandlimited so that
# it can be downsampled by 2 with no aliasing,
# i.e. X(e^jw)=0 for |w|>= pi/2
n = np.arange(length)
x = np.sin(np.pi/3*n)
# 'CT signal
ntrue = np.linspace(0,length,length*5)
xtrue = np.sin(np.pi/3*ntrue)
# remove some points
for idx in missing_points:
    x[idx] = -5
# interpolate to find the missing points
xintp = [elt for elt in x]
sinc = None
percent_error = []
for point in missing_points:
    # downsample then expand by 2
    # if the missing point had an even index, want to save all of the odd samples
    expanded = []
    if point%2==0:
        for w in range(len(x)):
            if w%2==0:
                expanded.append(0)
            else:
                expanded.append(x[w])            
    # if the missing point had an odd index, want to save all of the even samples
    elif point%2==1:
        for d in range(len(x)):
            if d%2==0:
                expanded.append(x[d])
            else:
                expanded.append(0)
    # windowed sinc filter
    #print '         forming windowed sinc filter...'
    L = 2
    n_neg = np.arange(-point,0)
    n_pos = np.arange(1,len(expanded)-point)
    sinc_neg = np.sin(np.pi*n_neg/L)/(np.pi*n_neg/L)
    sinc_pos = np.sin(np.pi*n_pos/L)/(np.pi*n_pos/L)
    # NB: normally would set sinc[0] equal to 1, but in this
    # case we don't want the -5 at that point to contribute to
    # the sum, so just set sinc[0] equal to 0
    sinc = np.concatenate((sinc_neg,np.array([0]),sinc_pos))
    # evaluate convolution at missing point
    missing = np.dot(expanded,sinc)
    #print '         missing point has value ',missing,' = ',np.sin(np.pi/3*point),' ?'
    xintp[point] = missing
    percent_error.append(abs(missing-math.sin(math.pi/3*point))/math.sin(math.pi/3*point)*100)

# plots
fig_xn = plt.figure()
# no interpolation
xn = fig_xn.add_subplot(211)
xn.stem(n[:20],x[:20],linestyle='none',marker='o')
xn.set_xlim(left=-1)
xn.set_ylim(bottom=-7,top=max(x[:20])*2)
# interpolated signal
xintpn = fig_xn.add_subplot(212)
xintpn.stem(n[:20],xintp[:20],linestyle='none',marker='o')
ntrue_20 = [elt for elt in ntrue if elt<20]
xtrue_20 = xtrue[:len(ntrue_20)]
xintpn.plot(ntrue_20,xtrue_20,linestyle='--')
xintpn.set_xlim(left=-1)
xintpn.set_ylim(bottom=-2,top=max(xintp[:20])*2)

# label with percent errors
fig_xn.text(0.90,0.93,str(percent_error[0])+'\n'+str(percent_error[1]),fontsize=10,color='red',va='baseline',ha='right',multialignment='left')
        
# set axis labels
x_axis_fontsize = 20
xintpn.set_xlabel('n',fontsize=x_axis_fontsize)

y_axis_fontsize = 20
xn.set_ylabel('x[n] \n (w/out interpolation)',fontsize=y_axis_fontsize-7)
xintpn.set_ylabel('x[n] \n (w/ interpolation)',fontsize=y_axis_fontsize-7)

# save figures
fig_xn.savefig('interpolation_sinusoid_lengthofsinc_'+str(length)+'.png')

plt.close('all')