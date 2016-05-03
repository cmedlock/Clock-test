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

path = '/Users/cmedlock/Documents/DSP_UROP/'
dirs = os.listdir(path)

pi = math.pi

n = np.arange(250)
omega = 2*math.pi/250.
axis_ratios = np.linspace(1,5,50)

# model order
p = 1
Eg_x = [] # energy in g[n] (i.e. linear prediction error) when predicting x[n]
Eg_y = [] # energy in g[n] (i.e. linear prediction error) when predicting y[n]

for ratio in axis_ratios:
    print 'ratio is now ',ratio
    x = np.cos(omega*n)
    y = ratio*np.sin(omega*n)
    
    # compensate for non-constant velocity
    
    N_orig = len(x)
    N_new = N_orig

    # calculate average distance between points
    dists = []
    for w in range(1,len(x)):
        dx,dy = x[w]-x[w-1],y[w]-y[w-1]
        dist = math.sqrt(dx**2+dy**2)
        dists.append(dist)
    dist_avg = mean(dists)
    dist_total = sum(dists)

    # if the points are already evenly-spaced, don't interpolate
    if np.var(np.array(dists))<10**-12:
        x_eqdist,y_eqdist = x,y
    else:
        # now want to get N_orig evenly-spaced points along the curve

        # generate a much longer array with 199 linearly-interpolated 
        # points between the actual data points
        x_interp,y_interp = [],[]
        for w in range(len(x)-1):
            x_interp.append(x[w])
            y_interp.append(y[w])
            dx,dy = x[w+1]-x[w],y[w+1]-y[w]
            dist = math.sqrt(dx**2+dy**2)
            n_segments = ceil(dist/dist_avg)*200
            for r in range(1,int(n_segments)):
                x_new = x[w]+r*dx/n_segments
                y_new = y[w]+r*dy/n_segments
                x_interp.append(x_new)
                y_interp.append(y_new)
        x_interp.append(x[-1])
        y_interp.append(y[-1])

        # start from the first point and find the ones that are 
        # approximately a distance dist_avg from each other
        x_eqdist,y_eqdist = [x_interp[0]],[y_interp[0]]
        idx = 0
        for k in range(N_new-1):
            dist_sofar = 0
            for j in range(idx,len(x_interp)-1):
                dx,dy = x_interp[j+1]-x_interp[j],y_interp[j+1]-y_interp[j]
                dist_sofar += math.sqrt(dx**2+dy**2)
                if abs(dist_sofar-dist_total/250.)<dist_total/(250.*100.):
                    idx = j+1
	            break
            x_eqdist.append(x_interp[idx])
            y_eqdist.append(y_interp[idx])

    # subtract mean values so there is no DC term adding an extra pole
    x_eqdist = [elt-mean(x_eqdist) for elt in x_eqdist]
    y_eqdist = [elt-mean(y_eqdist) for elt in y_eqdist]

    # downsample by a factor of 2
    x_downsampled,y_downsampled = [],[]
    for w in range(len(x_eqdist)):
        if w%2==0:
            x_downsampled.append(x_eqdist[w])
            y_downsampled.append(y_eqdist[w])
    #x_eqdist,y_eqdist = np.array(x_downsampled),np.array(y_downsampled)

    x_eqdist,y_eqdist = np.array(x_eqdist),np.array(y_eqdist)
    
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
fig.subplots_adjust(hspace=0.6,left=0.15)
ax = fig.add_subplot(111)
ax.clear()
ax.scatter(axis_ratios,Eg_x,color='blue',marker='o',s=100,alpha=0.5,label='x')
ax.scatter(axis_ratios,Eg_y,color='red',marker='^',s=100,alpha=0.5,label='y')
ax.set_xlabel('(major axis length)/(minor axis length)',fontsize=20)
ax.set_ylabel(r'$E_g/E_{total} \times 10^4$',fontsize=20)
ax.legend(loc='best')
fig.savefig(path+'circle_ellipse_lpc/compare_complex_Eg.png')
