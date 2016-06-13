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
p = 2
# save interesting quantities to compare between healthy and impaired patients
# format is [[('healthy',ak_1 val),('healthy',ak_1 val),('impaired',ak_1 val),...]
#            [ same thing for ak_2 ],...
#            [ same thing for ak_p]]
ak_x_coeffs,ak_y_coeffs = [],[]
Eg_x = [] # energy in g[n] (i.e. linear prediction error) when predicting x[n]
Eg_y = [] # energy in g[n] (i.e. linear prediction error) when predicting y[n]
for w in range(p):
    ak_x_coeffs.append([])
    ak_y_coeffs.append([])

axis_ratios = np.linspace(1,1.2,100)
for ratio in axis_ratios:

    x_eqdist = np.loadtxt(path+'norm velocity data/perfect_ellipse_x_eqdist_ratio_'+str(ratio)[:5]+'.txt')
    y_eqdist = np.loadtxt(path+'norm velocity data/perfect_ellipse_y_eqdist_ratio_'+str(ratio)[:5]+'.txt')
    
    # form all-pole model using Y-W eqns
    rxx,ryy = [],[]
    x_periodic,y_periodic = np.concatenate((x_eqdist,x_eqdist)),np.concatenate((y_eqdist,y_eqdist))
    for w in range(p+1):
        # circular autocorrelation method
        rxx.append(np.dot(x_eqdist,x_periodic[w:w+len(x_eqdist)]))
        ryy.append(np.dot(y_eqdist,y_periodic[w:w+len(y_eqdist)]))
    # calculate linear prediction coefficients
    D_x,D_y = np.array(rxx[1:p+1]),np.array(ryy[1:p+1])
    W_x,W_y = np.empty((p,p)),np.empty((p,p))
    ak_x,ak_y = np.empty((p)),np.empty((p))
    for row in range(p):
        for column in range(row,p):
            W_x[row][column] = rxx[column-row]
            W_x[column][row] = rxx[column-row]
            W_y[row][column] = ryy[column-row]
            W_y[column][row] = ryy[column-row]
    W_x_inv,W_y_inv = np.linalg.inv(W_x),np.linalg.inv(W_y)
    ak_x,ak_y = np.dot(W_x_inv,D_x),np.dot(W_y_inv,D_y)
    # linear prediction of x[n] or y[n]
    x_hat,y_hat = [],[]
    for w in range(len(x_eqdist)):
        prediction_x,prediction_y = 0,0
        for d in range(p):
            prediction_x += ak_x[d]*x_eqdist[w-(d+1)]
            prediction_y += ak_y[d]*y_eqdist[w-(d+1)]
        x_hat.append(prediction_x)
        y_hat.append(prediction_y)
    # linear prediction error of x[n] or y[n]
    g_x,g_y = x_eqdist-x_hat,y_eqdist-y_hat
    # percent energy in linear prediction error,
    energy_x = sum(g_x**2)/sum(np.array(x_eqdist)**2)*10**4
    energy_y = sum(g_y**2)/sum(np.array(y_eqdist)**2)*10**4
    
    # store the coefficients for comparison between the drawings of healthy
    # and impaired patients
    for m in range(p):
        ak_x_coeffs[m].append(ak_x[m])
        ak_y_coeffs[m].append(ak_y[m])
    Eg_x.append(energy_x)
    Eg_y.append(energy_y)

    # plot
    #plt.close('all')
    #np.set_printoptions(precision=2)
    #fig = plt.figure()
    ##ax = fig.add_subplot(111)
    ##ax.plot(x,y)
    ##fig.savefig(path+'circle_ellipse_lpc/rotated_ellipse_test_ratio_'+str(ratio)+'.png')
    #ax_lin,ax_lpe = fig.add_subplot(211),fig.add_subplot(212)
    ## predict x[n]
    #ax_lin.set_title('Linear Prediction of x[n]',fontsize=20)
    #ax_lin.plot(x_eqdist,color='blue',lw=2,label=r'$x_{eqdist}[n]$')
    #ax_lin.plot(x_hat,'k-.',lw=3,label=r'$\hat{x}[n]$')
    #ax_lin.legend(loc='best')
    #ax_lpe.plot(g_x,color='red',lw=2,label='$g[n]$')
    #ax_lpe.legend(loc='best')
    #fig.savefig(path+'circle_ellipse_lpc/lpc_x_ratio_'+str(ratio)+'.png')
    ## predict y[n]
    #ax_lin.clear()
    #ax_lpe.clear()
    #ax_lin.set_title('Linear Prediction of y[n]',fontsize=20)
    #ax_lin.plot(y_eqdist,color='blue',lw=2,label=r'$y_{eqdist}[n]$')
    #ax_lin.plot(y_hat,'k-.',lw=3,label=r'$\hat{y}[n]$')
    #ax_lin.legend(loc='best')
    #ax_lpe.plot(g_y,color='red',lw=2,label='$g[n]$')
    #ax_lpe.legend(loc='best')
    #fig.savefig(path+'circle_ellipse_lpc/lpc_y_ratio_'+str(ratio)+'.png')

# compare percent energy in linear prediction error
plt.close('all')
fig = plt.figure()
fig.subplots_adjust(hspace=0.6,left=0.2,bottom=0.15)
ax = fig.add_subplot(111)
ax.clear()

Eg_x_low_eccentricity,Eg_y_low_eccentricity = [],[]
Eg_x_high_eccentricity,Eg_y_high_eccentricity = [],[]
a1_x_low_eccentricity,a1_y_low_eccentricity = [],[]
a1_x_high_eccentricity,a1_y_high_eccentricity = [],[]
a2_x_low_eccentricity,a2_y_low_eccentricity = [],[]
a2_x_high_eccentricity,a2_y_high_eccentricity = [],[]
for w in range(len(axis_ratios)):
    if axis_ratios[w]<1.05:
        Eg_x_low_eccentricity.append(Eg_x[w])
        Eg_y_low_eccentricity.append(Eg_y[w])
        a1_x_low_eccentricity.append(ak_x_coeffs[0][w])
        a1_y_low_eccentricity.append(ak_y_coeffs[0][w])
        a2_x_low_eccentricity.append(ak_x_coeffs[1][w])
        a2_y_low_eccentricity.append(ak_y_coeffs[1][w])
    else:
        Eg_x_high_eccentricity.append(Eg_x[w])
        Eg_y_high_eccentricity.append(Eg_y[w])
        a1_x_high_eccentricity.append(ak_x_coeffs[0][w])
        a1_y_high_eccentricity.append(ak_y_coeffs[0][w])
        a2_x_high_eccentricity.append(ak_x_coeffs[1][w])
        a2_y_high_eccentricity.append(ak_y_coeffs[1][w])
"""
ax.scatter(Eg_x_low_eccentricity,Eg_y_low_eccentricity,color='green',marker='o',alpha=0.5,label='(major axis)/(minor axis) < 1.05')
ax.scatter(Eg_x_high_eccentricity,Eg_y_high_eccentricity,color='black',marker='o',alpha=0.5,label='(major axis)/(minor axis) > 1.05')
#ax.scatter(axis_ratios[:len(Eg_y_low_eccentricity)],Eg_x_low_eccentricity,color='green',marker='o',alpha=0.5,label='(major axis)/(minor axis) < 1.05')
#ax.scatter(axis_ratios[len(Eg_y_low_eccentricity):],Eg_x_high_eccentricity,color='black',marker='o',alpha=0.5,label='(major axis)/(minor axis) > 1.05')
ax.set_xlim(0.002,0.0035)
ax.set_ylim(0.8,1.05)
ax.set_xlabel(r'$E_g^x/E_{total} \times 10^4$',fontsize=20)
ax.set_ylabel(r'$E_g^y/E_{total} \times 10^4$',fontsize=20)
ax.legend(loc='best')
fig.savefig(path+'compare circles ellipses/compare_Eg_x_and_y.png')
"""
ax.clear()
ax.scatter(a1_x_low_eccentricity,a1_y_low_eccentricity,color='green',marker='o',alpha=0.5,label='(major axis)/(minor axis) < 1.05')
ax.scatter(a1_x_high_eccentricity,a1_y_high_eccentricity,color='black',marker='o',alpha=0.5,label='(major axis)/(minor axis) > 1.05')
ax.set_xlim(1.9991,1.99925)
ax.set_ylim(1.92,1.945)
y_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
ax.yaxis.set_major_formatter(y_formatter)
ax.xaxis.set_major_formatter(y_formatter)
ax.set_xlabel(r'$a_1$ for $x[n]$',fontsize=23)
ax.set_ylabel(r'$a_1$ for $y[n]$',fontsize=23)
j=1
for v1 in ax.get_xticklabels():
    v1.set_fontsize(12)
    if j==1:
        v1.set_visible(False)
    j += 1
j=1
for v2 in ax.get_yticklabels():
    v2.set_fontsize(17)
    #if j==1 or j%2==0:
    #    v1.set_visible(False)
    j += 1
ax.legend(loc='best')
fig.savefig(path+'compare circles ellipses/compare_a1_x_and_y.png')

ax.clear()
ax.scatter(a2_x_low_eccentricity,a2_y_low_eccentricity,color='green',marker='o',alpha=0.5,label='(major axis)/(minor axis) < 1.05')
ax.scatter(a2_x_high_eccentricity,a2_y_high_eccentricity,color='black',marker='o',alpha=0.5,label='(major axis)/(minor axis) > 1.05')
#ax.set_xlim(0.002,0.0035)
ax.set_ylim(-0.938,-0.923)
y_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
ax.yaxis.set_major_formatter(y_formatter)
ax.xaxis.set_major_formatter(y_formatter)
ax.set_xlabel(r'$a_2$ for $x[n]$',fontsize=23)
ax.set_ylabel(r'$a_2$ for $y[n]$',fontsize=23)
j=1
for v1 in ax.get_xticklabels():
    v1.set_fontsize(12)
    if j==1:
        v1.set_visible(False)
    j += 1
j=1
for v2 in ax.get_yticklabels():
    v2.set_fontsize(17)
    #if j==1 or j%2==0:
    #    v1.set_visible(False)
    j += 1
ax.legend(loc='best')
fig.savefig(path+'compare circles ellipses/compare_a2_x_and_y.png')
