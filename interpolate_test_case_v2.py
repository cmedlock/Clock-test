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
# this code only works for even lengths
lengths = [20,40]
missing_point = 5

for length in lengths:
    print '***length is now ',length,'***'
    # original signal, appropriately bandlimited so that
    # it can be downsampled by 2 with no aliasing,
    # i.e. X(e^jw)=0 for |w|>= pi/2
    n = np.arange(missing_point-length/2,missing_point+length/2)
    x = np.sin(np.pi/6*n)
    # 'CT signal
    ntrue = np.linspace(missing_point-length/2,missing_point+length/2,length*5)
    xtrue = np.sin(np.pi/6*ntrue)
    # remove a point
    idx_of_missing_point = -1
    for j in range(len(n)):
        idx = n[j]
        if idx==missing_point:
            x[j] = -5
            print 'idx of missing point is ',j
            idx_of_missing_point = j
    # interpolate to find the missing points
    xintp = [elt for elt in x]
    sinc = []
    percent_error = -1
    
    # downsample then expand by 2
    # if the missing point had an even index, want to save all of the odd samples
    expanded = []
    if missing_point%2==0:
        for w in range(len(x)):
            if w%2==0:
                expanded.append(0)
            else:
                expanded.append(x[w])            
    # if the missing point had an odd index, want to save all of the even samples
    elif missing_point%2==1:
        for d in range(len(x)):
            if d%2==0:
                expanded.append(x[d])
            else:
                expanded.append(0)
    # windowed sinc filter
    #print '         forming windowed sinc filter...'
    L = 2
    n_neg = n[:idx_of_missing_point]
    n_pos = n[idx_of_missing_point+1:]
    n_neg,n_pos = n_neg-n_neg[-1]-1,n_pos-n_pos[0]+1
    sinc_neg = np.sin(np.pi*n_neg/L)/(np.pi*n_neg/L)
    sinc_pos = np.sin(np.pi*n_pos/L)/(np.pi*n_pos/L)    
    # NB: normally would set sinc[0] equal to 1, but in this
    # case we don't want the -5 at that point to contribute to
    # the sum, so just set sinc[0] equal to 0
    sinc = np.concatenate((sinc_neg,np.array([0]),sinc_pos))
    print '   ',x[idx_of_missing_point],sinc[idx_of_missing_point]
    # evaluate convolution at missing point
    missing = np.dot(expanded,sinc)
    print '   missing point has value ',missing,' = ',math.sin(math.pi/6*missing_point),' ?'
    xintp[idx_of_missing_point] = missing
    #print x[idx_of_missing_point],xintp[idx_of_missing_point]
    #error = abs(missing-math.sin(math.pi/6*(point-40)))/math.sin(math.pi/6*(point-length/2))*100
    #print 'percent error for point ',point,' is ',error
    #percent_error.append(error)

    # plots
    fig_xn = plt.figure()
    # no interpolation
    xn = fig_xn.add_subplot(311)
    xn.stem(x,linestyle='none',marker='o')
    xn.set_xlim(left=-1)
    xn.set_ylim(bottom=-7,top=max(x)*2)
    # downsampled-then-expanded signal and truncated sinc
    sincn = fig_xn.add_subplot(312)
    sincn.stem(sinc,linestyle='none',marker='o')
    sincn.set_xlim(left=-1)
    # interpolated signal
    xintpn = fig_xn.add_subplot(313)
    xintpn.stem(xintp,linestyle='none',marker='o')
    #xintpn.plot(ntrue,xtrue,linestyle='--')
    xintpn.set_xlim(left=-1)
    xintpn.set_ylim(bottom=-2,top=max(xintp)*2)
    
    # label with percent errors
    #fig_xn.text(0.90,0.93,str(percent_error[0])+'\n'+str(percent_error[1]),fontsize=10,color='red',va='baseline',ha='right',multialignment='left')
            
    # set axis labels
    x_axis_fontsize = 20
    xintpn.set_xlabel('n',fontsize=x_axis_fontsize)
    
    y_axis_fontsize = 20
    xn.set_ylabel('x[n] \n (w/out interpolation)',fontsize=y_axis_fontsize-7)
    sincn.set_ylabel('sinc[n]',fontsize=y_axis_fontsize-7)
    xintpn.set_ylabel('x[n] \n (w/ interpolation)',fontsize=y_axis_fontsize-7)
    
    # save figures
    fig_xn.savefig('interpolation_sinusoid_lengthofsinc_'+str(length)+'.png')
        
plt.close('all')
