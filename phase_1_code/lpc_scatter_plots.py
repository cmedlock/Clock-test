# -*- coding: utf-8 -*-
# 1.find LPC cepstral coefficients as described in signature analysis paper
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

path = '/Users/cmedlock/Documents/DSP_UROP/DataCatherine/'
dirs = os.listdir(path)

pi = math.pi

if not os.path.exists(path+'figs_raw'):
    os.makedirs(path+'figs_raw')

def make_scatter_plots(clock_type):
    # linear prediction coefficients
    for m in range(2):
        pvals_x_healthy = np.loadtxt('pvals_x_healthy_a'+str(m+1)+'_'+clock_type+'.txt')
        pvals_x_impaired = np.loadtxt('pvals_x_impaired_a'+str(m+1)+'_'+clock_type+'.txt')
    
        ak_x_healthy = np.loadtxt('a'+str(m+1)+'_x_healthy_'+clock_type+'.txt')
        ak_x_impaired = np.loadtxt('a'+str(m+1)+'_x_impaired_'+clock_type+'.txt')
    
        #plt.close('all')
        #fig = plt.figure()
        #ax = fig.add_subplot(111)
        #ax.clear()
        #ax.scatter(pvals_x_healthy,ak_x_healthy,color='red',marker='o',s=200,alpha=0.5,label='Healthy')
        #ax.scatter(pvals_x_impaired,ak_x_impaired,color='blue',marker='o',s=200,alpha=0.5,label='Impaired')
        #ax.legend(loc='lower center',frameon=False)
        #ax.set_xlabel('p',fontsize=20)
        #ax.set_ylabel('a_'+str(m+1)+' (x)',fontsize=20)
        #fig.savefig(path+'compare_healthy_impaired/compare_a'+str(m+1)+'_x_'+clock_type+'.png')
        if m==0:
            range_min,range_max = 0.9,1.1
        elif m==1:
            range_min,range_max = -0.1,0.1
        ct.make_hist(ak_x_healthy,ak_x_impaired,
                20,range_min,range_max,'a_'+str(m+1)+' for x[n]','a'+str(m+1)+'_x_'+clock_type,path)

        pvals_y_healthy = np.loadtxt('pvals_y_healthy_a'+str(m+1)+'_'+clock_type+'.txt')
        pvals_y_impaired = np.loadtxt('pvals_y_impaired_a'+str(m+1)+'_'+clock_type+'.txt')
        
        ak_y_healthy = np.loadtxt('a'+str(m+1)+'_y_healthy_'+clock_type+'.txt')
        ak_y_impaired = np.loadtxt('a'+str(m+1)+'_y_impaired_'+clock_type+'.txt')
    
        #ax.clear()
        #ax.scatter(pvals_y_healthy,ak_y_healthy,color='red',marker='o',s=200,alpha=0.5,label='Healthy')
        #ax.scatter(pvals_y_impaired,ak_y_impaired,color='blue',marker='o',s=200,alpha=0.5,label='Impaired')
        #ax.legend(loc='lower center',frameon=False)
        #ax.set_xlabel('p',fontsize=20)
        #ax.set_ylabel('a_'+str(m+1)+' (y)',fontsize=20)
        #fig.savefig(path+'compare_healthy_impaired/compare_a'+str(m+1)+'_y_'+clock_type+'.png')
        if m==0:
            range_min,range_max = 0,2
        elif m==1:
            range_min,range_max = -1,0
        ct.make_hist(ak_y_healthy,ak_y_impaired,
                20,range_min,range_max,'a_'+str(m+1)+' for y[n]','a'+str(m+1)+'_y_'+clock_type,path)

    # energy in linear prediction error compared to energy of signal
    # modeling x[n] using x[n]
    pvals_healthy = np.loadtxt('pvals_healthy_'+clock_type+'.txt')
    pvals_impaired = np.loadtxt('pvals_impaired_'+clock_type+'.txt')

    Eg_x_healthy = np.loadtxt('Eg_x_healthy_'+clock_type+'.txt')
    Eg_x_impaired = np.loadtxt('Eg_x_impaired_'+clock_type+'.txt')

    #ax.clear()
    #ax.scatter(pvals_healthy,Eg_x_healthy,color='red',marker='o',s=200,alpha=0.5,label='Healthy')
    #ax.scatter(pvals_impaired,Eg_x_impaired,color='blue',marker='o',s=200,alpha=0.5,label='Impaired')
    #ax.legend(loc='lower center',frameon=False)
    #ax.set_xlabel('p',fontsize=20)
    #ax.set_ylabel(r'$E_g/E_{total} \times 10^2 = $',fontsize=20)
    #ax.set_title('Relative Energy in Linear Prediction Error for x[n]',fontsize=30)
    #fig.savefig(path+'compare_healthy_impaired/compare_LPE_x_'+clock_type+'.png')
    ct.make_hist(Eg_x_healthy,Eg_x_impaired,
                20,0.05,0.15,'Relative Energy in Linear Prediction Error for x[n]','LPE_x_'+clock_type,path)

    # modeling y[n] using y[n]
    Eg_y_healthy = np.loadtxt('Eg_y_healthy_'+clock_type+'.txt')
    Eg_y_impaired = np.loadtxt('Eg_y_impaired_'+clock_type+'.txt')

    #ax.clear()
    #ax.scatter(pvals_healthy,Eg_y_healthy,color='red',marker='o',s=200,alpha=0.5,label='Healthy')
    #ax.scatter(pvals_impaired,Eg_y_impaired,color='blue',marker='o',s=200,alpha=0.5,label='Impaired')
    #ax.legend(loc='lower center',frameon=False)
    #ax.set_xlabel('p',fontsize=20)
    #ax.set_ylabel(r'$E_g/E_{total} \times 10^2 = $',fontsize=20)
    #ax.set_title('Relative Energy in Linear Prediction Error for y[n]',fontsize=30)
    #fig.savefig(path+'compare_healthy_impaired/compare_LPE_y_'+clock_type+'.png')
    ct.make_hist(Eg_y_healthy,Eg_y_impaired,
                20,0,0.5,'Relative Energy in Linear Prediction Error for y[n]','LPE_y_'+clock_type,path)

# copy or command clock?
clock_type = 'COMMAND'

# model order
pvals = [2]
# save interesting quantities to compare between healthy and impaired patients
# format is [[('healthy',ak_1 val),('healthy',ak_1 val),('impaired',ak_1 val),...]
#            [ same thing for ak_2 ],...
#            [ same thing for ak_p]]
ak_x_coeffs,ak_y_coeffs = [],[]
Eg_x_x,Eg_x_y = [],[] # energy in g[n] (i.e. linear prediction error) when predicting x[n]
Eg_y_x,Eg_y_y = [],[] # energy in g[n] (i.e. linear prediction error) when predicting y[n]
for w in range(3):
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
        continue
        
    # subtract mean values so there is no DC term adding an extra pole
    x_eqdist = np.array(x_eqdist)-mean(x_eqdist)
    y_eqdist = np.array(y_eqdist)-mean(y_eqdist)

    # form all-pole model using Y-W eqns
    for p in pvals:
        # covariance circular autocorrelation method
        rxx,ryy = [],[]
        x_periodic,y_periodic = np.concatenate((x_eqdist,x_eqdist)),np.concatenate((y_eqdist,y_eqdist))
        for w in range(p+1):
            if not cov:
                # circular autocorrelation method
                rxx.append(np.dot(x_eqdist,x_periodic[w:w+len(x_eqdist)]))
                ryy.append(np.dot(y_eqdist,y_periodic[w:w+len(y_eqdist)]))
            else:
                # covariance method
                rxx.append(np.dot(x_eqdist[p-w:],x_periodic[p-w:len(x_eqdist)]))
                ryy.append(np.dot(y_eqdist[p-w:],y_periodic[p-w:len(y_eqdist)]))
        # calculate linear prediction coefficients
        D_x,D_y = np.array(rxx[1:p+1]),np.array(ryy[1:p+1])
        W_x,W_y = np.empty((p,p)),np.empty((p,p))
        ak_x,ak_y = np.empty((p)),np.empty((p))
        ck_x,ck_y = [],[]
        for row in range(p):
            for column in range(row,p):
                W_x[row][column] = rxx[column-row]
                W_x[column][row] = rxx[column-row]
                W_y[row][column] = ryy[column-row]
                W_y[column][row] = ryy[column-row]
        # LPC spectrum
        W_x_inv,W_y_inv = np.linalg.inv(W_x),np.linalg.inv(W_y)
        ak_x,ak_y = np.dot(W_x_inv,D_x),np.dot(W_y_inv,D_y)
        # impulse response of linear prediction filter
        a_x,a_y = [1],[1]
        for w in range(p):
            a_x.append(-float(ak_x[w]))
            a_y.append(-float(ak_y[w]))

        # linear prediction error (LPE):
        g_x_x = np.empty((len(x_eqdist)+p)) # predict x[n] with x[n] model
        g_x_y = np.empty((len(x_eqdist)+p)) # predict x[n] with y[n] model
        g_y_x = np.empty((len(y_eqdist)+p)) # predict y[n] with x[n] model
        g_y_y = np.empty((len(y_eqdist)+p)) # predict y[n] with y[n] model
        for w in range(len(x_eqdist)+p):
            lpe_x_x,lpe_x_y,lpe_y_x,lpe_y_y = 0,0,0,0
            #if w==15:
                #print 'w = ',w
            for d in range(len(a_x)):
                if w-d<len(x_eqdist)+p:
                    #if w==15:
                        #print '***   ',w-d,d,' ---> ',x_periodic[w-d],a_x[d]
                    lpe_x_x += x_periodic[w-d]*a_x[d]
                    lpe_x_y += x_periodic[w-d]*a_y[d]
                    lpe_y_x += x_periodic[w-d]*a_x[d]
                    lpe_y_y += y_periodic[w-d]*a_y[d]
            g_x_x[w] = lpe_x_x
            g_x_y[w] = lpe_x_y
            g_y_x[w] = lpe_y_x
            g_y_y[w] = lpe_y_y
        g_x_x,g_x_y,g_y_x,g_y_y = np.array(g_x_x),np.array(g_x_y),np.array(g_y_x),np.array(g_y_y)
        # total energy in linear prediction error,
        # normalized by total energy in original signal
        energy_x_x = sum(g_x_x[:-p]**2)/sum(np.array(x_eqdist)**2)*10**2
        energy_x_y = sum(g_x_y[:-p]**2)/sum(np.array(x_eqdist)**2)*10**2
        energy_y_x = sum(g_y_x[:-p]**2)/sum(np.array(y_eqdist)**2)*10**2
        energy_y_y = sum(g_y_y[:-p]**2)/sum(np.array(y_eqdist)**2)*10**2
        # linear prediction of x[n] or y[n]
        xhat_x_x = -g_x_x[:-p]+x_eqdist
        xhat_x_y = -g_x_y[:-p]+x_eqdist
        yhat_y_x = -g_y_x[:-p]+y_eqdist
        yhat_y_y = -g_y_y[:-p]+y_eqdist
    
        # store the coefficients for comparison between the drawings of healthy
        # and impaired patients
        for m in range(p):
            ak_x_coeffs[m].append((ftype,(p,ak_x[m])))
            ak_y_coeffs[m].append((ftype,(p,ak_y[m])))
        Eg_x_x.append((ftype,(p,energy_x_x)))
        Eg_x_y.append((ftype,(p,energy_x_y)))
        Eg_y_x.append((ftype,(p,energy_y_x)))
        Eg_y_y.append((ftype,(p,energy_y_y)))

        ## plot
        #plt.close('all')
        #np.set_printoptions(precision=2)
        #fig = plt.figure()
        #ax_lin,ax_lpe = fig.add_subplot(211),fig.add_subplot(212)
        #fig.text(0.99, 0.96,fname[:len(fname)-4],fontsize=10,color='red',va='baseline',ha='right',multialignment='left')
        #if ftype=='healthy':
        #    fig.text(0.32, 0.955, 'HEALTHY ('+clock_type+')',fontsize=15,color='black',va='baseline',ha='right',multialignment='left')
        #elif ftype=='impaired':
        #    fig.text(0.32, 0.955, 'IMPAIRED ('+clock_type+')',fontsize=15,color='black',va='baseline',ha='right',multialignment='left')
        ## predict x[n]
        #ax_lin.set_title('Linear Prediction of x[n]',fontsize=20)
        #ax_lin.plot(x_eqdist,color='blue',lw=2,label=r'$x_{eqdist}[n]$')
        #ax_lin.plot(xhat_x_x,'k-.',lw=3,label=r'$\hat{x}[n]$')
        #ax_lin.legend(loc='best')
        #ax_lpe.plot(g_x_x[:-p],color='red',lw=2,label='$g[n]$')
        #ax_lpe.legend(loc='best')
        #ax_lpe.text(10,1.5,r'$E_g/E_{total} \times 10^2 = $'+str(np.round_(energy_x_x,3)),fontsize=15)
        #ax_lpe.text(10,2.3,'p = '+str(p),fontsize=15)
        #ax_lin.set_ylim(bottom=-70,top=70)
        #ax_lpe.set_ylim(bottom=-3,top=3)
        #fig.savefig(path+'figs_raw/'+fname[:len(fname)-4]+'/lpc_x_p'+str(p)+'_'+clock_type+'_'+fname[:len(fname)-4]+'.png')
        ## predict y[n]
        #ax_lin.clear()
        #ax_lpe.clear()
        #ax_lin.set_title('Linear Prediction of y[n]',fontsize=20)
        #ax_lin.plot(y_eqdist,color='blue',lw=2,label=r'$y_{eqdist}[n]$')
        #ax_lin.plot(yhat_y_y,'k-.',lw=3,label=r'$\hat{y}[n]$')
        #ax_lin.legend(loc='best')
        #ax_lpe.plot(g_y_y[:-p],color='red',lw=2,label='$g[n]$')
        #ax_lpe.legend(loc='best')
        #ax_lpe.text(10,2.3,'p = '+str(p),fontsize=15)
        #ax_lpe.text(10,1.5,r'$E_g/E_{total} \times 10^2 = $'+str(np.round_(energy_y_y,3)),fontsize=15)
        #ax_lin.set_ylim(bottom=-70,top=70)
        #ax_lpe.set_ylim(bottom=-3,top=3)
        #fig.savefig(path+'figs_raw/'+fname[:len(fname)-4]+'/lpc_y_p'+str(p)+'_'+clock_type+'_'+fname[:len(fname)-4]+'.png')

# compare linear prediction coefficients and total
# energy in linear prediction error between the drawings of healthy
# and impaired patients

for m in range(2):
    pvals_x_healthy = [elt[1][0] for elt in ak_x_coeffs[m] if elt[0]=='healthy']
    pvals_x_impaired = [elt[1][0] for elt in ak_x_coeffs[m] if elt[0]=='impaired']
    np.savetxt('pvals_x_healthy_a'+str(m+1)+'_'+clock_type+'.txt',pvals_x_healthy)
    np.savetxt('pvals_x_impaired_a'+str(m+1)+'_'+clock_type+'.txt',pvals_x_impaired)
    
    ak_x_healthy = [elt[1][1] for elt in ak_x_coeffs[m] if elt[0]=='healthy']
    ak_x_impaired = [elt[1][1] for elt in ak_x_coeffs[m] if elt[0]=='impaired']
    np.savetxt('a'+str(m+1)+'_x_healthy_'+clock_type+'.txt',ak_x_healthy)
    np.savetxt('a'+str(m+1)+'_x_impaired_'+clock_type+'.txt',ak_x_impaired)

    pvals_y_healthy = [elt[1][0] for elt in ak_y_coeffs[m] if elt[0]=='healthy']
    pvals_y_impaired = [elt[1][0] for elt in ak_y_coeffs[m] if elt[0]=='impaired']
    np.savetxt('pvals_y_healthy_a'+str(m+1)+'_'+clock_type+'.txt',pvals_y_healthy)
    np.savetxt('pvals_y_impaired_a'+str(m+1)+'_'+clock_type+'.txt',pvals_y_impaired)
   
    ak_y_healthy = [elt[1][1] for elt in ak_y_coeffs[m] if elt[0]=='healthy']
    ak_y_impaired = [elt[1][1] for elt in ak_y_coeffs[m] if elt[0]=='impaired']
    np.savetxt('a'+str(m+1)+'_y_healthy_'+clock_type+'.txt',ak_y_healthy)
    np.savetxt('a'+str(m+1)+'_y_impaired_'+clock_type+'.txt',ak_y_impaired)

pvals_healthy = [elt[1][0] for elt in Eg_x_x if elt[0]=='healthy']
pvals_impaired = [elt[1][0] for elt in Eg_x_x if elt[0]=='impaired']
np.savetxt('pvals_healthy_'+clock_type+'.txt',pvals_healthy)
np.savetxt('pvals_impaired_'+clock_type+'.txt',pvals_impaired)

Eg_x_healthy = [elt[1][1] for elt in Eg_x_x if elt[0]=='healthy']
Eg_x_impaired = [elt[1][1] for elt in Eg_x_x if elt[0]=='impaired']
np.savetxt('Eg_x_healthy_'+clock_type+'.txt',Eg_x_healthy)
np.savetxt('Eg_x_impaired_'+clock_type+'.txt',Eg_x_impaired)

Eg_y_healthy = [elt[1][1] for elt in Eg_y_y if elt[0]=='healthy']
Eg_y_impaired = [elt[1][1] for elt in Eg_y_y if elt[0]=='impaired']
np.savetxt('Eg_y_healthy_'+clock_type+'.txt',Eg_y_healthy)
np.savetxt('Eg_y_impaired_'+clock_type+'.txt',Eg_y_impaired)

make_scatter_plots(clock_type)
