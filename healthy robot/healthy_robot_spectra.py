# for sanity checks on spectral and energy properties

import math
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import numpy.fft
import os
from pylab import *

import clock_test as ct
ct = reload(ct)
    
path = '/Users/cmedlock/Documents/DSP_UROP/all_data/'
dirs = os.listdir(path)

if not os.path.exists(path+'figs_raw'):
    os.makedirs(path+'figs_raw')

for fname in dirs[:3]:
    
    if 'Scored' not in fname:
        continue
    
    if not os.path.exists(path+'figs_raw/'+fname[:len(fname)-4]):
        os.makedirs(path+'figs_raw/'+fname[:len(fname)-4])

    # get coordinates and timestamps
    n = np.arange(250)
    x = np.cos(2*np.pi/250*n)
    y = np.sin(2*np.pi/250*n)
        
    # plot spectrum of x[t] and y[t]
    # subtract mean values (zm = zero mean)
    x_zm,y_zm = x-mean(x),y-mean(y)
    # dft's
    dft_size = len(n)
    dftx,dfty = np.fft.fft(x_zm,n=dft_size),np.fft.fft(y_zm,n=dft_size)
    print max(np.abs(dftx)),np.abs(dftx[1])
    print max(np.abs(dfty)),np.abs(dfty[1])
    k = np.arange(dft_size)
    freq = 2*np.pi/dft_size*k
    
    # k_near_pi is the smallest k value for which w_k = 2*pi*k/N is
    # greater than or equal to pi
    k_near_pi = 0
    if dft_size%2==0:
        k_near_pi = dft_size/2+1
    else:
        k_near_pi = math.ceil(dft_size/2)

    k_centered = np.linspace(-dft_size/2,dft_size/2,dft_size)
    dftx_centered = np.concatenate((dftx[k_near_pi:],dftx[:k_near_pi]))
    dfty_centered = np.concatenate((dfty[k_near_pi:],dfty[:k_near_pi]))

    # percent energy in peak
    Ex,Ey = np.abs(dftx_centered)**2,np.abs(dfty_centered)**2
    Ex_total,Ey_total = sum(Ex),sum(Ey)
    Ex_peak,Ey_peak = 2*max(Ex)/Ex_total,2*max(Ey)/Ey_total
    print 'Ex_peak = ',Ex_peak,' and Ey_peak = ',Ey_peak
    
    # percent energy within 1 std. deviation of the center of the energy distribution
    Ex_var,Ey_var = 0,0
    for w in range(dft_size):
        Ex_var += k_centered[w]**2*Ex[w]/Ex_total
        Ey_var += k_centered[w]**2*Ey[w]/Ey_total
    Ex_std,Ey_std = math.sqrt(Ex_var),math.sqrt(Ey_var)
    print 'Ex_std = ',Ex_std,' and Ey_std = ',Ey_std
    
    Ex_central,Ey_central = 0,0
    for d in range(dft_size):
        if abs((dft_size-1)/2-d)<=Ex_std:
            Ex_central += Ex[d]/Ex_total
        if abs((dft_size-1)/2-d)<=Ey_std:
            Ey_central += Ey[d]/Ey_total
    print 'Ex_central = ',Ex_central,' and Ey_central = ',Ey_central

    # plot
    ct.plot_xyt_other(x,n,'x','n',np.abs(dftx)[:k_near_pi],k[:k_near_pi],'|X[k]|','k',True,'dftx','healthy_robot',path)
    ct.plot_xyt_other(y,n,'y','n',np.abs(dfty)[:k_near_pi],k[:k_near_pi],'|Y[k]|','k',True,'dfty','healthy_robot',path)

    plt.close('all')
    
print 'Done'
