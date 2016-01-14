# -*- coding: utf-8 -*-
import math
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import numpy.fft
import os
from pylab import *

path = '/Users/cmedlock/Documents/DSP_UROP/all_data/'

np.set_printoptions(precision=2)

if not os.path.exists(path+'figs_raw'):
    os.makedirs(path+'figs_raw')

# length of truncated sinc (longer = more accurate interpolation)
# this code only works for even lengths
lengths = [10,20,60,100,200,300,400,500,600,700,800,900,1000]
#lengths = [10]
missing_point = 1
percent_errors = []

for length in lengths:
    check = 0
    for k in range(missing_point-length/2,missing_point+length/2):
        if k%2==0:
            check += math.sin(math.pi/6*k)*math.sin(math.pi/2*(k-1))/(math.pi/2*(k-1))
            #print k,'   ',math.sin(math.pi/6*k),'   ',math.sin(math.pi/2*(k-1))/(math.pi/2*(k-1))

    print '***length is now ',length,'***'
    # original signal, appropriately bandlimited so that
    # it can be downsampled by 2 with no aliasing,
    # i.e. X(e^jw)=0 for |w|>= pi/2
    n = np.arange(missing_point-length/2,missing_point+length/2)
    #print 'n is ',n
    x = np.sin(np.pi/6*n)
    # 'CT signal
    ntrue = np.linspace(missing_point-length/2,missing_point+length/2,length*5)
    xtrue = np.sin(np.pi/6*ntrue)
    # remove a point
    idx_of_missing_point = -1
    for j in range(len(n)):
        idx = n[j]
        if idx==missing_point:
            #print 'removing point ',j,',',x[j]
            x[j] = -5
            idx_of_missing_point = j
    # interpolate to find the missing points
    xintp = [elt for elt in x]
    sinc = []
    
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
        for d in range(len(n)):
            nval = n[d]
            if nval%2==0:
                expanded.append(x[d])
            else:
                expanded.append(0)
    # check
    for w in range(len(n)):
        nval = n[w]
        if nval%2==0:
            pass
            #print nval,expanded[w],math.sin(math.pi/6*nval)

    # windowed sinc filter
    #print '         forming windowed sinc filter...'
    L = 2
    n_neg = n[:idx_of_missing_point]
    n_pos = n[idx_of_missing_point+1:]
    n_neg,n_pos = n_neg-n_neg[-1]-1,n_pos-n_pos[0]+1
    #print 'n_neg and n_pos are ',n_neg,n_pos
    sinc_neg = np.sin(np.pi/L*n_neg)/(np.pi/L*n_neg)
    sinc_pos = np.sin(np.pi/L*n_pos)/(np.pi/L*n_pos)    
    # NB: normally would set sinc[0] equal to 1, but in this
    # case we don't want the -5 at that point to contribute to
    # the sum, so just set sinc[0] equal to 0
    sinc = np.concatenate((sinc_neg,np.array([0]),sinc_pos))
    # check
    for w in range(len(n_neg)):
        nval = n_neg[w]
        pass
        #print nval,': ',sinc_neg[w],'   ',math.sin(math.pi/2*nval)/(math.pi/2*nval)
    # evaluate convolution at missing point
    missing = np.dot(expanded,sinc)
    #print 'missing point has value ',missing,' = ',check,' = ',math.sin(math.pi/6*missing_point),' ?'
    xintp[idx_of_missing_point] = missing
    error = abs(missing-math.sin(math.pi/6*missing_point))/math.sin(math.pi/6*missing_point)*100
    print '      percent error for point ',missing_point,' is ',error
    percent_errors.append(error)

    # plots
    if length<100:
        fig_xn = plt.figure()
        # no interpolation
        xn = fig_xn.add_subplot(311)
        xn.stem(n,x,linestyle='none',marker='o')
        xn.set_xlim(left=n[0]-0.5,right=n[-1]+0.5)
        xn.set_ylim(bottom=-7,top=max(x)*2)
        # downsampled-then-expanded signal and truncated sinc
        sincn = fig_xn.add_subplot(312)
        sincn.stem(n,sinc,linestyle='none',marker='o')
        sincn.set_xlim(left=n[0]-0.5,right=n[-1]+0.5)
        # interpolated signal
        xintpn = fig_xn.add_subplot(313)
        xintpn.stem(n,xintp,linestyle='none',marker='o')
        #xintpn.plot(ntrue,xtrue,linestyle='--')
        xintpn.set_xlim(left=n[0]-0.5,right=n[-1]+0.5)
        xintpn.set_ylim(bottom=-2,top=max(xintp)*2)
    
        # label with percent error
        fig_xn.text(0.90,0.93,'percent error = '+str(error),fontsize=15,color='red',va='baseline',ha='right',multialignment='left')
            
        # set axis labels
        x_axis_fontsize = 20
        xintpn.set_xlabel('n',fontsize=x_axis_fontsize)
    
        y_axis_fontsize = 20
        xn.set_ylabel('x[n] \n (w/out interpolation)',fontsize=y_axis_fontsize-7)
        sincn.set_ylabel('sinc[n]',fontsize=y_axis_fontsize-7)
        xintpn.set_ylabel('x[n] \n (w/ interpolation)',fontsize=y_axis_fontsize-7)
    
        # save figures
        fig_xn.savefig('interpolation_sinusoid_lengthofsinc_'+str(length)+'.png')

# percent error plot
fig_percent_error = plt.figure()
# no interpolation
percent_err = fig_percent_error.add_subplot(111)
percent_err.plot(lengths,percent_errors,marker='o')
percent_err.set_xlim(right=max(lengths)*1.05)
percent_err.set_ylim(top=max(percent_errors)*1.2)
    
# set axis labels
x_axis_fontsize = 20
percent_err.set_xlabel('length of truncated sinc',fontsize=x_axis_fontsize)
    
y_axis_fontsize = 20
percent_err.set_ylabel(r'$\%$'+' error',fontsize=y_axis_fontsize)
    
# save figures
fig_percent_error.savefig('interpolation_percent_error.png')
        
plt.close('all')