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

def sinc(omega_c,n,length_of_sinc):
    if n==0:
        return 1
    elif abs(n)>(length_of_sinc-1)/2:
        return 0
    else:
        return math.sin(omega_c*n)/(omega_c*n)

path = '/Users/cmedlock/Documents/DSP_UROP/all_data/'
dirs = os.listdir(path)

pi = math.pi

if not os.path.exists(path+'figs_raw'):
    os.makedirs(path+'figs_raw')

# copy or command clock?
clock_type = 'COPY'

# model order
pvals = [2,3]
# covariance method?
cov = True
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
    Eg_x_x.append([])
    Eg_x_y.append([])
    Eg_y_x.append([])
    Eg_y_y.append([])

for fname in dirs:
    if 'Scored' not in fname:
        continue
    print 'reading file ',fname,'...'
    
    ftype = ''
    if 'YDU' in fname:
        ftype = 'healthy'
    elif 'CIN' in fname:
        ftype = 'impaired'
    else:
        print 'not a valid file name'

    f = open(path+fname)
    data = f.readlines()

    if not os.path.exists(path+'figs_raw/'+fname[:len(fname)-4]):
        os.makedirs(path+'figs_raw/'+fname[:len(fname)-4])

    x,y,t = [],[],[]

    # read in data
    record,found_clock = False,False
    for w in range(len(data)):
        line = data[w]
        # found clock?
        if found_clock==False:
            if clock_type in line:
                found_clock = True
            continue
        # start recording?
        elif found_clock==True and 'CLOCKFACE' in line:
            record = True
            continue
        # stop recording?
        elif record==True:
            if 'symbol label' in line and len(x)>0:
                record = False
                break
            elif 'point' not in line:
                continue
        # done?
        elif found_clock==False and record==False and len(x)>0:
            break
        # other?
        else:
            continue
        line = line.split('"')
        
        # position and time
        xcoord = double(line[3])
        ycoord = double(line[1])
        timestamp = double(line[7])
        x.append(xcoord)
        y.append(ycoord)
        t.append(timestamp)
    
    f.close()

    # change x so that if the whole clock is drawn,
    # it is oriented correctly
    x = max(x)+10-x
    
    # normalize the size of the word such that all coordinates
    # are linearly positioned between 0 and 127
    x = 127*(x-min(x))/(max(x)-min(x))
    y = 127*(y-min(y))/(max(y)-min(y))

    # compensate for non-constant velocity
    
    N_orig = len(x)
    N_new = 250

    # calculate average distance between points
    dists = []
    for w in range(1,len(x)):
        dx,dy = x[w]-x[w-1],y[w]-y[w-1]
        dist = math.sqrt(dx**2+dy**2)
        dists.append(dist)
    dist_avg = mean(dists)
    dist_total = sum(dists)
    #print 'average distance between points is ',dist_avg_copy
    #print 'total distance is ',sum(dists_copy)

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
    for k in range(N_new):
        dist_sofar = 0
        for j in range(idx,len(x_interp)-1):
            dx,dy = x_interp[j+1]-x_interp[j],y_interp[j+1]-y_interp[j]
            dist_sofar += math.sqrt(dx**2+dy**2)
            #if abs(dist_sofar-dist_avg)<dist_avg/100.:
            if abs(dist_sofar-dist_total/250.)<dist_total/(250.*100.):
                idx = j+1
	        break
        x_eqdist.append(x_interp[idx])
        y_eqdist.append(y_interp[idx])

    # subtract mean values so there is no DC term adding an 'extra' pole
    x_eqdist = [elt-mean(x_eqdist) for elt in x_eqdist]
    y_eqdist = [elt-mean(y_eqdist) for elt in y_eqdist]
    
    # form all-pole model using Y-W eqns
    for p in pvals:
        # circular autocorrelation
        # or covariance
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
            g_x_x[w] = lpe_x_x if abs(lpe_x_x)>10**-15 else 0
            g_x_y[w] = lpe_x_y if abs(lpe_x_y)>10**-15 else 0
            g_y_x[w] = lpe_y_x if abs(lpe_y_x)>10**-15 else 0
            g_y_y[w] = lpe_y_y if abs(lpe_y_y)>10**-15 else 0
        g_x_x,g_x_y,g_y_x,g_y_y = np.array(g_x_x),np.array(g_x_y),np.array(g_y_x),np.array(g_y_y)
        # total energy in linear prediction error,
        # normalized by total energy in original signal
        # multiply by 10**3 so that the numbers are
        # readable
        energy_x_x = sum(g_x_x[:-p]**2)/sum(np.array(x_eqdist)**2)*10**3
        energy_x_y = sum(g_x_y[:-p]**2)/sum(np.array(x_eqdist)**2)*10**3
        energy_y_x = sum(g_y_x[:-p]**2)/sum(np.array(y_eqdist)**2)*10**3
        energy_y_y = sum(g_y_y[:-p]**2)/sum(np.array(y_eqdist)**2)*10**3
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
            Eg_x_x[m].append((ftype,(p,energy_x_x)))
            Eg_x_y[m].append((ftype,(p,energy_x_y)))
            Eg_y_x[m].append((ftype,(p,energy_y_x)))
            Eg_y_y[m].append((ftype,(p,energy_y_y)))
        
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
        # predict x[n] with x[n] model
        ax_lin.set_title('Model x[n] using x[n]',fontsize=20)
        ax_lin.plot(x_eqdist,color='blue',lw=2,label=r'$x_{eqdist}[n]$')
        ax_lin.plot(xhat_x_x,'k-.',lw=3,label=r'$\hat{x}[n]$')
        ax_lin.legend(loc='best')
        ax_lpe.plot(g_x_x[:-p],color='red',lw=2,label='$g[n]$')
        ax_lpe.legend(loc='best')
        ax_lpe.text(10,1.5,r'$E_g/E_{total} \times 10^3 = $'+str(np.round_(energy_x_x,3)),fontsize=15)
        ax_lpe.text(10,2.3,'p = '+str(p),fontsize=15)
        ax_lin.set_ylim(bottom=-70,top=70)
        ax_lpe.set_ylim(bottom=-3,top=3)
        fig.savefig(path+'figs_raw/'+fname[:len(fname)-4]+'/lpc_x_x_p'+str(p)+'_cov'+str(cov)+'_'+clock_type+'_'+fname[:len(fname)-4]+'.png')
        # predict x[n] with y[n] model
        ax_lin.clear()
        ax_lpe.clear()
        ax_lin.set_title('Model x[n] using y[n]',fontsize=20)
        ax_lin.plot(x_eqdist,color='blue',lw=2,label=r'$x_{eqdist}[n]$')
        ax_lin.plot(xhat_x_y,'k-.',lw=3,label=r'$\hat{x}[n]$')
        ax_lin.legend(loc='best')
        ax_lpe.plot(g_x_y[:-p],color='red',lw=2,label='$g[n]$')
        ax_lpe.legend(loc='best')
        ax_lpe.text(10,2.3,'p = '+str(p),fontsize=15)
        ax_lpe.text(10,1.5,r'$E_g/E_{total} \times 10^3 = $'+str(np.round_(energy_x_y,3)),fontsize=15)
        ax_lin.set_ylim(bottom=-70,top=70)
        ax_lpe.set_ylim(bottom=-3,top=3)
        fig.savefig(path+'figs_raw/'+fname[:len(fname)-4]+'/lpc_x_y_p'+str(p)+'_cov'+str(cov)+'_'+clock_type+'_'+fname[:len(fname)-4]+'.png')
        # predict y[n] with x[n] model
        ax_lin.clear()
        ax_lpe.clear()
        ax_lin.set_title('Model y[n] using x[n]',fontsize=20)
        ax_lin.plot(y_eqdist,color='blue',lw=2,label=r'$y_{eqdist}[n]$')
        ax_lin.plot(yhat_y_x,'k-.',lw=3,label=r'$\hat{y}[n]$')
        ax_lin.legend(loc='best')
        ax_lpe.plot(g_y_x[:-p],color='red',lw=2,label='$g[n]$')
        ax_lpe.legend(loc='best')
        ax_lpe.text(10,2.3,'p = '+str(p),fontsize=15)
        ax_lpe.text(10,1.5,r'$E_g/E_{total} \times 10^3 = $'+str(np.round_(energy_y_x,3)),fontsize=15)
        ax_lin.set_ylim(bottom=-70,top=70)
        ax_lpe.set_ylim(bottom=-3,top=3)
        fig.savefig(path+'figs_raw/'+fname[:len(fname)-4]+'/lpc_y_x_p'+str(p)+'_cov'+str(cov)+'_'+clock_type+'_'+fname[:len(fname)-4]+'.png')
        # predict y[n] with y[n] model
        ax_lin.clear()
        ax_lpe.clear()
        ax_lin.set_title('Model y[n] using y[n]',fontsize=20)
        ax_lin.plot(y_eqdist,color='blue',lw=2,label=r'$y_{eqdist}[n]$')
        ax_lin.plot(yhat_y_y,'k-.',lw=3,label=r'$\hat{y}[n]$')
        ax_lin.legend(loc='best')
        ax_lpe.plot(g_y_y[:-p],color='red',lw=2,label='$g[n]$')
        ax_lpe.legend(loc='best')
        ax_lpe.text(10,2.3,'p = '+str(p),fontsize=15)
        ax_lpe.text(10,1.5,r'$E_g/E_{total} \times 10^3 = $'+str(np.round_(energy_y_y,3)),fontsize=15)
        ax_lin.set_ylim(bottom=-70,top=70)
        ax_lpe.set_ylim(bottom=-3,top=3)
        fig.savefig(path+'figs_raw/'+fname[:len(fname)-4]+'/lpc_y_y_p'+str(p)+'_cov'+str(cov)+'_'+clock_type+'_'+fname[:len(fname)-4]+'.png')

# compare linear prediction coefficients and total
# energy in linear prediction error between the drawings of healthy
# and impaired patients
plt.close('all')
fig = plt.figure()
ax = fig.add_subplot(111)
# linear prediction coefficients
for m in range(2):
    pvals_x_healthy = [elt[1][0] for elt in ak_x_coeffs[m] if elt[0]=='healthy']
    pvals_x_impaired = [elt[1][0] for elt in ak_x_coeffs[m] if elt[0]=='impaired']
    
    ak_x_healthy = [elt[1][1] for elt in ak_x_coeffs[m] if elt[0]=='healthy']
    ak_x_impaired = [elt[1][1] for elt in ak_x_coeffs[m] if elt[0]=='impaired']
    
    ax.clear()
    ax.plot(pvals_x_healthy,ak_x_healthy,color='red',marker='*')
    ax.set_xlabel('p',fontsize=20)
    ax.set_ylabel('a_'+str(m+1)+' (x)',fontsize=20)
    fig.savefig(path+'compare_healthy_impaired/compare_a'+str(m+1)+'_p'+str(p)+'_cov'+str(cov)+'_'+clock_type+'.png')
