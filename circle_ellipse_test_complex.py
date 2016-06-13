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

# path to figures for each individual file
if not os.path.exists(path+'compare circles ellipses'):
    os.makedirs(path+'compare circles ellipses')
    
pi = math.pi

# model order
p = 1
Eg_x = [] # energy in g[n] (i.e. linear prediction error) when predicting x[n]
Eg_y = [] # energy in g[n] (i.e. linear prediction error) when predicting y[n]

axis_ratios = np.linspace(1,1.2,100)
for ratio in axis_ratios:

    x_eqdist = np.loadtxt(path+'norm velocity data/perfect_ellipse_x_eqdist_ratio_'+str(ratio)[:5]+'.txt')
    y_eqdist = np.loadtxt(path+'norm velocity data/perfect_ellipse_y_eqdist_ratio_'+str(ratio)[:5]+'.txt')
    
    z_eqdist = x_eqdist+1j*y_eqdist
    z_eqdist = z_eqdist+0j

    # form all-pole model using Y-W eqns
    rzz = []
    z_periodic = np.concatenate((z_eqdist,z_eqdist))
    for w in range(p+1):
        # circular autocorrelation method
        rzz.append(np.vdot(z_eqdist,z_periodic[w:w+len(z_eqdist)]))
    # calculate linear prediction coefficients
    W_z = rzz[0]
    D_z = rzz[1]
    ak_z = D_z/W_z
    # linear prediction
    z_hat = np.array([0j]*len(z_eqdist),dtype=np.complex)
    z_hat = z_hat+0j
    for w in range(len(z_eqdist)):
        prediction_z = ak_z*z_eqdist[w-1]
        z_hat[w] = prediction_z
    x_hat = np.real(z_hat)
    y_hat = np.imag(z_hat)
    # linear prediction error
    g_x = x_eqdist-x_hat
    g_y = y_eqdist-y_hat
    # percent energy in linear prediction error,
    energy_x = sum(g_x**2)/sum(np.array(x_eqdist)**2)*10**4
    energy_y = sum(g_y**2)/sum(np.array(y_eqdist)**2)*10**4
    Eg_x.append(energy_x)
    Eg_y.append(energy_y)
    
    # plot
    plt.close('all')
    #np.set_printoptions(precision=2)
    #fig = plt.figure()
    #ax_lin,ax_lpe = fig.add_subplot(211),fig.add_subplot(212)
    ## predict x[n]
    #ax_lin.set_title('Linear Prediction of x[n]',fontsize=20)
    #ax_lin.plot(x_eqdist,color='blue',lw=2,label=r'$x_{eqdist}[n]$')
    #ax_lin.plot(x_hat,'k-.',lw=3,label=r'$\hat{x}[n]$')
    #ax_lin.legend(loc='best')
    #ax_lpe.plot(g_x,color='red',lw=2,label='$g[n]$')
    #ax_lpe.legend(loc='best')
    #fig.savefig(path+'circle_ellipse_lpc/complex_lpc_x_ratio_'+str(ratio)+'.png')
    ## predict y[n]
    #ax_lin.clear()
    #ax_lpe.clear()
    #ax_lin.set_title('Linear Prediction of y[n]',fontsize=20)
    #ax_lin.plot(y_eqdist,color='blue',lw=2,label=r'$y_{eqdist}[n]$')
    #ax_lin.plot(y_hat,'k-.',lw=3,label=r'$\hat{y}[n]$')
    #ax_lin.legend(loc='best')
    #ax_lpe.plot(g_y,color='red',lw=2,label='$g[n]$')
    #ax_lpe.legend(loc='best')
    #fig.savefig(path+'circle_ellipse_lpc/complex_lpc_y_ratio_'+str(ratio)+'.png')

# compare percent energy in linear prediction error
plt.close('all')
fig = plt.figure()
fig.subplots_adjust(hspace=0.6,left=0.15,bottom=0.15)
ax = fig.add_subplot(111)
ax.clear()

Eg_x_low_eccentricity,Eg_y_low_eccentricity = [],[]
Eg_x_high_eccentricity,Eg_y_high_eccentricity = [],[]
for w in range(len(axis_ratios)):
    if axis_ratios[w]<1.05:
        Eg_x_low_eccentricity.append(Eg_x[w])
        Eg_y_low_eccentricity.append(Eg_y[w])
    else:
        Eg_x_high_eccentricity.append(Eg_x[w])
        Eg_y_high_eccentricity.append(Eg_y[w])

ax.scatter(Eg_x_low_eccentricity,Eg_y_low_eccentricity,color='green',marker='o',alpha=0.5,label='(major axis)/(minor axis) < 1.05')
ax.scatter(Eg_x_high_eccentricity,Eg_y_high_eccentricity,color='black',marker='o',alpha=0.5,label='(major axis)/(minor axis) > 1.05')
#ax.set_xlim(0,2)
#ax.set_ylim(0.035,0.05)
ax.set_xlabel(r'$E_g^x/E_{total} \times 10^4$',fontsize=20)
ax.set_ylabel(r'$E_g^y/E_{total} \times 10^4$',fontsize=20)
for v1 in ax.get_xticklabels():
    v1.set_fontsize(17)
for v2 in ax.get_yticklabels():
    v2.set_fontsize(17)
ax.legend(loc='lower right')
fig.savefig(path+'compare circles ellipses/compare_complex_Eg_x_and_y.png')
"""
ax.clear()
Eg_x = np.concatenate((Eg_x_low_eccentricity,Eg_x_high_eccentricity))
Eg_y = np.concatenate((Eg_y_low_eccentricity,Eg_y_high_eccentricity))
ax.scatter(axis_ratios,Eg_x,color='blue',marker='o',alpha=0.5,label=r'$E_g^x$')
ax.scatter(axis_ratios,Eg_y,color='red',marker='o',alpha=0.5,label=r'$E_g^y$')
ax.set_xlabel('(major axis)/(minor axis)',fontsize=20)
ax.set_ylabel(r'$E_g/E_{total} \times 10^4$',fontsize=20)
ax.legend(loc='upper left')
fig.savefig(path+'circle_ellipse_lpc/compare_complex_Eg.png')
"""