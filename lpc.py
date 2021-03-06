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

path = '/Users/cmedlock/Documents/DSP_UROP/data_for_report/'
dirs = os.listdir(path)

pi = math.pi

# copy or command clock?
clock_type = 'copy'

# path to linear prediction figures
if not os.path.exists(path+'lpc'):
    os.makedirs(path+'lpc')
if not os.path.exists(path+'lpc/'+clock_type+'_healthy'):
    os.makedirs(path+'lpc/'+clock_type+'_healthy')
if not os.path.exists(path+'lpc/'+clock_type+'_impaired'):
    os.makedirs(path+'lpc/'+clock_type+'_impaired')
if not os.path.exists(path+'compare_healthy_impaired'):
    os.makedirs(path+'compare_healthy_impaired')

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
    energy_x = sum(g_x**2)/sum(np.array(x_eqdist)**2)*10**2
    energy_y = sum(g_y**2)/sum(np.array(y_eqdist)**2)*10**2
    
    # store the coefficients for comparison between the drawings of healthy
    # and impaired patients
    for m in range(p):
        ak_x_coeffs[m].append((ftype,ak_x[m]))
        ak_y_coeffs[m].append((ftype,ak_y[m]))
    Eg_x.append((ftype,energy_x))
    Eg_y.append((ftype,energy_y))

    # plot
    plt.close('all')
    np.set_printoptions(precision=2)
    fig = plt.figure()
    ax_lin,ax_lpe = fig.add_subplot(211),fig.add_subplot(212)
    #fig.text(0.99, 0.96,fname[:len(fname)-4],fontsize=10,color='red',va='baseline',ha='right',multialignment='left')
    #if ftype=='healthy':
    #    fig.text(0.32, 0.955, 'HEALTHY ('+clock_type+')',fontsize=15,color='black',va='baseline',ha='right',multialignment='left')
    #elif ftype=='impaired':
    #    fig.text(0.32, 0.955, 'IMPAIRED ('+clock_type+')',fontsize=15,color='black',va='baseline',ha='right',multialignment='left')
    # predict x[n]
    #ax_lin.set_title('Linear Prediction of x[n]',fontsize=20)
    ax_lin.plot(x_eqdist,color='blue',lw=2,label=r'$x_{eqdist}[n]$')
    ax_lin.plot(x_hat,'k-.',lw=3,label=r'$\hat{x}[n]$')
    #ax_lin.legend(loc='best')
    ax_lpe.plot(g_x,color='red',lw=2,label='$g[n]$')
    #ax_lpe.legend(loc='best')
    # unify axis limits
    y_min_lin,y_max_lin = min(min(x_eqdist),min(y_eqdist)),max(max(x_eqdist),max(y_eqdist))
    y_min_lpe,y_max_lpe = min(min(g_x),min(g_y)),max(max(g_x),max(g_y))
    ax_lin.set_ylim(y_min_lin-5,y_max_lin+5)
    ax_lpe.set_ylim(y_min_lpe-0.1,y_max_lpe+0.1)
    j = 1
    for v1 in ax_lin.get_yticklabels():
        v1.set_fontsize(21)
        #if j%2==0:
        #    v1.set_visible(False)
        j += 1
    j = 1
    for v2 in ax_lpe.get_yticklabels():
        v2.set_fontsize(21)
        if j%2==0:
            v2.set_visible(False)
        j += 1
    j = 1
    for v1 in ax_lin.get_xticklabels():
        v1.set_fontsize(21)
        if j%2==0:
            v1.set_visible(False)
        j += 1
    j = 1
    for v2 in ax_lpe.get_xticklabels():
        v2.set_fontsize(21)
        if j%2==0:
            v2.set_visible(False)
        j += 1
    fig.savefig(path+'lpc/'+clock_type+'_'+ftype+'/lpc_x_'+clock_type+'_'+fname[:len(fname)-4]+'.png')
    fig.savefig(path+'figs_raw/'+fname[:len(fname)-4]+'/lpc_x_'+clock_type+'_'+fname[:len(fname)-4]+'.png')
    # predict y[n]
    ax_lin.clear()
    ax_lpe.clear()
    #ax_lin.set_title('Linear Prediction of y[n]',fontsize=20)
    ax_lin.plot(y_eqdist,color='blue',lw=2,label=r'$y_{eqdist}[n]$')
    ax_lin.plot(y_hat,'k-.',lw=3,label=r'$\hat{y}[n]$')
    #ax_lin.legend(loc='best')
    ax_lpe.plot(g_y,color='red',lw=2,label='$g[n]$')
    #ax_lpe.legend(loc='best')
    # unify axis limits
    ax_lin.set_ylim(y_min_lin-5,y_max_lin+5)
    ax_lpe.set_ylim(y_min_lpe-0.1,y_max_lpe+0.1)
    j = 1
    for v1 in ax_lin.get_yticklabels():
        v1.set_fontsize(21)
        #if j%2==0:
        #    v1.set_visible(False)
        j += 1
    j = 1
    for v2 in ax_lpe.get_yticklabels():
        v2.set_fontsize(21)
        #if j%2==0:
        #    v2.set_visible(False)
        j += 1
    j = 1
    for v1 in ax_lin.get_xticklabels():
        v1.set_fontsize(21)
        if j%2==0:
            v1.set_visible(False)
        j += 1
    j = 1
    for v2 in ax_lpe.get_xticklabels():
        v2.set_fontsize(21)
        if j%2==0:
            v2.set_visible(False)
        j += 1
    fig.savefig(path+'lpc/'+clock_type+'_'+ftype+'/lpc_y_'+clock_type+'_'+fname[:len(fname)-4]+'.png')
    fig.savefig(path+'figs_raw/'+fname[:len(fname)-4]+'/lpc_y_'+clock_type+'_'+fname[:len(fname)-4]+'.png')

# compare linear prediction coefficients and total
# energy in linear prediction error between the drawings of healthy and impaired patients
Eg_x_all = [elt[1] for elt in Eg_x]
mean_x,std_x = mean(Eg_x_all),std(Eg_x_all)
Eg_y_all = [elt[1] for elt in Eg_y]
mean_y,std_y = mean(Eg_y_all),std(Eg_y_all)
#ct.make_hist([elt[1] for elt in Eg_x if elt[0]=='healthy'],
#             [elt[1] for elt in Eg_x if elt[0]=='impaired'],
#             15,mean_x-std_x,mean_x+std_x,'Relative Energy in Linear Prediction Error for x[n]','Eg_x_'+clock_type,path)
#ct.make_hist([elt[1] for elt in Eg_y if elt[0]=='healthy'],
#             [elt[1] for elt in Eg_y if elt[0]=='impaired'],
#             15,mean_y-std_y,mean_y+std_y,'Relative Energy in Linear Prediction Error for y[n]','Eg_y_'+clock_type,path)
# in case the histograms don't come out right
#np.savetxt(path+'lpc/Eg_x_healthy_'+clock_type+'.txt',[elt[1] for elt in Eg_x if elt[0]=='healthy'])
#np.savetxt(path+'lpc/Eg_x_impaired_'+clock_type+'.txt',[elt[1] for elt in Eg_x if elt[0]=='impaired'])
#np.savetxt(path+'lpc/Eg_y_healthy_'+clock_type+'.txt',[elt[1] for elt in Eg_y if elt[0]=='healthy'])
#np.savetxt(path+'lpc/Eg_y_impaired_'+clock_type+'.txt',[elt[1] for elt in Eg_y if elt[0]=='impaired'])
np.savetxt(path+'lpc/a1_x_healthy_'+clock_type+'.txt',[elt[1] for elt in ak_x_coeffs[0] if elt[0]=='healthy'])
np.savetxt(path+'lpc/a1_x_impaired_'+clock_type+'.txt',[elt[1] for elt in ak_x_coeffs[0] if elt[0]=='impaired'])
np.savetxt(path+'lpc/a1_y_healthy_'+clock_type+'.txt',[elt[1] for elt in ak_y_coeffs[0] if elt[0]=='healthy'])
np.savetxt(path+'lpc/a1_y_impaired_'+clock_type+'.txt',[elt[1] for elt in ak_y_coeffs[0] if elt[0]=='impaired'])
np.savetxt(path+'lpc/a2_x_healthy_'+clock_type+'.txt',[elt[1] for elt in ak_x_coeffs[1] if elt[0]=='healthy'])
np.savetxt(path+'lpc/a2_x_impaired_'+clock_type+'.txt',[elt[1] for elt in ak_x_coeffs[1] if elt[0]=='impaired'])
np.savetxt(path+'lpc/a2_y_healthy_'+clock_type+'.txt',[elt[1] for elt in ak_y_coeffs[1] if elt[0]=='healthy'])
np.savetxt(path+'lpc/a2_y_impaired_'+clock_type+'.txt',[elt[1] for elt in ak_y_coeffs[1] if elt[0]=='impaired'])
