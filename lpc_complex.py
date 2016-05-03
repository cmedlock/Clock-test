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

path = '/Users/cmedlock/Documents/DSP_UROP/all_data/'
dirs = os.listdir(path)

pi = math.pi

# copy or command clock?
clock_type = 'command'

# path to linear prediction figures
if not os.path.exists(path+'lpc_complex'):
    os.makedirs(path+'lpc_complex')
if not os.path.exists(path+'lpc_complex/'+clock_type+'_healthy'):
    os.makedirs(path+'lpc_complex/'+clock_type+'_healthy')
if not os.path.exists(path+'lpc_complex/'+clock_type+'_impaired'):
    os.makedirs(path+'lpc_complex/'+clock_type+'_impaired')
if not os.path.exists(path+'compare_healthy_impaired'):
    os.makedirs(path+'compare_healthy_impaired')

# model order
p = 1
Eg_x = [] # energy in g[n] (i.e. linear prediction error) when predicting x[n]
Eg_y = [] # energy in g[n] (i.e. linear prediction error) when predicting y[n]

for fname in dirs:
    if 'CIN' not in fname and 'YDU' not in fname:
        continue
    print 'reading file ',fname,'...'
    
    ftype = ''
    if 'YDU' in fname:
        ftype = 'healthy'
    elif 'CIN' in fname:
        ftype = 'impaired'
    else:
        print 'not a valid file name'

    x_eqdist = np.loadtxt(path+'norm_velocity_data/'+fname[:len(fname)-4]+'_x_eqdist_'+clock_type+'.txt')
    y_eqdist = np.loadtxt(path+'norm_velocity_data/'+fname[:len(fname)-4]+'_y_eqdist_'+clock_type+'.txt')
    if len(x_eqdist)==0:
        print fname
        continue
    
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
    energy_x = sum(g_x**2)/sum(np.array(x_eqdist)**2)*10**2
    energy_y = sum(g_y**2)/sum(np.array(y_eqdist)**2)*10**2
    Eg_x.append((ftype,energy_x))
    Eg_y.append((ftype,energy_y))

    # plot
    plt.close('all')
    np.set_printoptions(precision=2)
    fig = plt.figure()
    ax_lin,ax_lpe = fig.add_subplot(211),fig.add_subplot(212)
    fig.text(0.99, 0.96,fname[:len(fname)-4],fontsize=10,color='red',va='baseline',ha='right',multialignment='left')
    if ftype=='healthy':
        fig.text(0.32, 0.955, 'HEALTHY ('+clock_type+')',fontsize=15,color='black',va='baseline',ha='right',multialignment='left')
    elif ftype=='impaired':
        fig.text(0.32, 0.955, 'IMPAIRED ('+clock_type+')',fontsize=15,color='black',va='baseline',ha='right',multialignment='left')
    # predict x[n]
    ax_lin.set_title('Linear Prediction of x[n]',fontsize=20)
    ax_lin.plot(x_eqdist,color='blue',lw=2,label=r'$x_{eqdist}[n]$')
    ax_lin.plot(x_hat,'k-.',lw=3,label=r'$\hat{x}[n]$')
    ax_lin.legend(loc='best')
    ax_lpe.plot(g_x,color='red',lw=2,label='$g[n]$')
    ax_lpe.legend(loc='best')
    fig.savefig(path+'lpc_complex/'+clock_type+'_'+ftype+'/lpc_complex_x_'+clock_type+'_'+fname[:len(fname)-4]+'.png')
    fig.savefig(path+'figs_raw/'+fname[:len(fname)-4]+'/lpc_complex_x_'+clock_type+'_'+fname[:len(fname)-4]+'.png')
    # predict y[n]
    ax_lin.clear()
    ax_lpe.clear()
    ax_lin.set_title('Linear Prediction of y[n]',fontsize=20)
    ax_lin.plot(y_eqdist,color='blue',lw=2,label=r'$y_{eqdist}[n]$')
    ax_lin.plot(y_hat,'k-.',lw=3,label=r'$\hat{y}[n]$')
    ax_lin.legend(loc='best')
    ax_lpe.plot(g_y,color='red',lw=2,label='$g[n]$')
    ax_lpe.legend(loc='best')
    fig.savefig(path+'lpc_complex/'+clock_type+'_'+ftype+'/lpc_complex_y_'+clock_type+'_'+fname[:len(fname)-4]+'.png')
    fig.savefig(path+'figs_raw/'+fname[:len(fname)-4]+'/lpc_complex_y_'+clock_type+'_'+fname[:len(fname)-4]+'.png')

# compare percent energy in linear prediction error between the drawings of healthy and impaired patients
Eg_x_all = [elt[1] for elt in Eg_x]
mean_x,std_x = mean(Eg_x_all),std(Eg_x_all)
Eg_y_all = [elt[1] for elt in Eg_y]
mean_y,std_y = mean(Eg_y_all),std(Eg_y_all)
ct.make_hist([elt[1] for elt in Eg_x if elt[0]=='healthy'],
             [elt[1] for elt in Eg_x if elt[0]=='impaired'],
             15,mean_x-std_x,mean_x+std_x,'Relative Energy in Linear Prediction Error for x[n]','complex_Eg_x_'+clock_type,path)
ct.make_hist([elt[1] for elt in Eg_y if elt[0]=='healthy'],
             [elt[1] for elt in Eg_y if elt[0]=='impaired'],
             15,mean_y-std_y,mean_y+std_y,'Relative Energy in Linear Prediction Error for y[n]','complex_Eg_y_'+clock_type,path)
# in case the histograms don't come out right
np.savetxt(path+'lpc_complex/Eg_x_healthy_'+clock_type+'.txt',[elt[1] for elt in Eg_x if elt[0]=='healthy'])
np.savetxt(path+'lpc_complex/Eg_x_impaired_'+clock_type+'.txt',[elt[1] for elt in Eg_x if elt[0]=='impaired'])
np.savetxt(path+'lpc_complex/Eg_y_healthy_'+clock_type+'.txt',[elt[1] for elt in Eg_y if elt[0]=='healthy'])
np.savetxt(path+'lpc_complex/Eg_y_impaired_'+clock_type+'.txt',[elt[1] for elt in Eg_y if elt[0]=='impaired'])
