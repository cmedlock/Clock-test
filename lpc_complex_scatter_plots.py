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

def make_hists(clock_type):
    error_x_healthy = np.loadtxt('complex_lpc_error_x_healthy_'+clock_type+'.txt')
    error_x_impaired = np.loadtxt('complex_lpc_error_x_impaired_'+clock_type+'.txt')
    error_y_healthy = np.loadtxt('complex_lpc_error_y_healthy_'+clock_type+'.txt')
    error_y_impaired = np.loadtxt('complex_lpc_error_y_impaired_'+clock_type+'.txt')
    # energy in linear prediction error compared to energy of signal
    ct.make_hist(error_x_healthy,error_x_impaired,
                20,'Relative Energy in Linear Prediciton Error (Complex)','complex_Eg_x_'+clock_type,path)
    ct.make_hist(error_y_healthy,error_y_impaired,
                20,'Relative Energy in Linear Prediciton Error (Complex)','complex_Eg_y_'+clock_type,path)

path = '/Users/cmedlock/Documents/DSP_UROP/DataCatherine/'
dirs = os.listdir(path)

pi = math.pi

if not os.path.exists(path+'figs_raw'):
    os.makedirs(path+'figs_raw')

# copy or command clock?
clock_type = 'COPY'
#
## model order
#pvals = [1]
## save interesting quantities to compare between healthy and impaired patients
## format is [[('healthy',ak_1 val),('healthy',ak_1 val),('impaired',ak_1 val),...]
##            [ same thing for ak_2 ],...
##            [ same thing for ak_p]]
#Eg_z_x,Eg_z_y = [],[] # energy in g[n] (i.e. linear prediction error)
#
#for fname in dirs:
#    if 'CIN' not in fname and 'YDU' not in fname:
#        continue
#    print 'reading file ',fname,'...'
#    
#    ftype = ''
#    if 'YDU' in fname:
#        ftype = 'healthy'
#    elif 'CIN' in fname:
#        ftype = 'impaired'
#    else:
#        print 'not a valid file name'
#
#    f = open(path+fname)
#    data = f.readlines()
#
#    if not os.path.exists(path+'figs_raw/'+fname[:len(fname)-4]):
#        os.makedirs(path+'figs_raw/'+fname[:len(fname)-4])
#
#    x,y,t = [],[],[]
#
#    # read in data
#    record,found_clock = False,False
#    for w in range(len(data)):
#        line = data[w]
#        # found clock?
#        if found_clock==False:
#            if clock_type in line:
#                found_clock = True
#            continue
#        # start recording?
#        elif found_clock==True and 'CLOCKFACE' in line:
#            record = True
#            continue
#        # stop recording?
#        elif record==True:
#            if 'symbol label' in line and len(x)>0:
#                record = False
#                break
#            elif 'point' not in line:
#                continue
#        # done?
#        elif found_clock==False and record==False and len(x)>0:
#            break
#        # other?
#        else:
#            continue
#        line = line.split('"')
#        
#        # position and time
#        xcoord = double(line[3])
#        ycoord = double(line[1])
#        timestamp = double(line[7])
#        x.append(xcoord)
#        y.append(ycoord)
#        t.append(timestamp)
#    
#    f.close()
#    if len(x)==0:
#        continue
#
#    # change x so that if the whole clock is drawn,
#    # it is oriented correctly
#    #x = max(x)+10-x
#    
#    # normalize the size of the word such that all coordinates
#    # are linearly positioned between 0 and 127
#    #x = 127*(x-min(x))/(max(x)-min(x))
#    #y = 127*(y-min(y))/(max(y)-min(y))
#
#    # compensate for non-constant velocity
#    
#    N_orig = len(x)
#    N_new = 250
#
#    # calculate average distance between points
#    dists = []
#    for w in range(1,len(x)):
#        dx,dy = x[w]-x[w-1],y[w]-y[w-1]
#        dist = math.sqrt(dx**2+dy**2)
#        dists.append(dist)
#    dist_avg = mean(dists)
#    dist_total = sum(dists)
#    #print 'average distance between points is ',dist_avg_copy
#    #print 'total distance is ',sum(dists_copy)
#
#    # if the points are already evenly-spaced, don't interpolate
#    if np.var(np.array(dists))<10**-12:
#        x_eqdist,y_eqdist = x,y
#    else:
#        # now want to get N_orig evenly-spaced points along the curve
#
#        # generate a much longer array with 199 linearly-interpolated 
#        # points between the actual data points
#        x_interp,y_interp = [],[]
#        for w in range(len(x)-1):
#            x_interp.append(x[w])
#            y_interp.append(y[w])
#            dx,dy = x[w+1]-x[w],y[w+1]-y[w]
#            dist = math.sqrt(dx**2+dy**2)
#            n_segments = ceil(dist/dist_avg)*200
#            for r in range(1,int(n_segments)):
#                x_new = x[w]+r*dx/n_segments
#                y_new = y[w]+r*dy/n_segments
#                x_interp.append(x_new)
#                y_interp.append(y_new)
#        x_interp.append(x[-1])
#        y_interp.append(y[-1])
#
#        # start from the first point and find the ones that are 
#        # approximately a distance dist_avg from each other
#        x_eqdist,y_eqdist = [x_interp[0]],[y_interp[0]]
#        idx = 0
#        for k in range(N_new-1):
#            dist_sofar = 0
#            for j in range(idx,len(x_interp)-1):
#                dx,dy = x_interp[j+1]-x_interp[j],y_interp[j+1]-y_interp[j]
#                dist_sofar += math.sqrt(dx**2+dy**2)
#                #if abs(dist_sofar-dist_avg)<dist_avg/100.:
#                if abs(dist_sofar-dist_total/250.)<dist_total/(250.*100.):
#                    idx = j+1
#	            break
#                x_eqdist.append(x_interp[idx])
#                y_eqdist.append(y_interp[idx])
#
#    # subtract mean values so there is no DC term adding an extra pole
#    x_eqdist = np.array(x_eqdist)-mean(x_eqdist)
#    y_eqdist = np.array(y_eqdist)-mean(y_eqdist)
#    
#    # form complex signal
#    z_eqdist = x_eqdist+1j*y_eqdist
#    z_eqdist = z_eqdist+0j
#
#    # form all-pole model using Y-W eqns
#    for p in pvals:
#        # circular autocorrelation method
#        rzz = []
#        z_periodic = np.concatenate((z_eqdist,z_eqdist))
#        for w in range(p+1):
#            # circular autocorrelation method
#            rzz.append(np.vdot(z_eqdist,z_periodic[w:w+len(z_eqdist)]))
#        # calculate linear prediction coefficients
#        if p==1:
#            W_z = rzz[0]
#            D_z = rzz[1]
#        elif p==2:
#            W_z = np.array([[0j,0j],[0j,0j]])
#            D_z = np.array(rzz[1:p+1])
#        ak_z = np.array([0j]*p)
#        for row in range(p):
#            for column in range(row,p):
#                if p==2:
#                    W_z[row][column] = rzz[column-row]
#                    W_z[column][row] = rzz[column-row]
#        # LPC spectrum
#        if p==1:
#            ak_z = np.array([D_z/W_z])
#        elif p==2:
#            W_z_inv = np.linalg.inv(W_z)
#            ak_z = np.dot(W_z_inv,D_z)
#        # impulse response of linear prediction filter
#        a_z = [1+0j]
#        for w in range(p):
#            a_z.append(-ak_z[w])
#
#        # linear prediction error (LPE):
#        g_z = np.empty((len(z_eqdist)+p))
#        g_z = g_z+0j
#        for w in range(len(z_eqdist)+p):
#            lpe_z = 0+0j
#            for d in range(len(a_z)):
#                if w-d<len(z_eqdist)+p:
#                    lpe_z += z_periodic[w-d]*a_z[d]
#            g_z[w] = lpe_z
#        g_z = np.array(g_z)
#        # linear prediction of x[n] or y[n]
#        zhat = -g_z[:-p]+z_eqdist
#        zhat_x = np.real(zhat)
#        zhat_y = np.imag(zhat)
#        gz_x = x_eqdist-zhat_x
#        gz_y = y_eqdist-zhat_y
#        # total energy in linear prediction error,
#        # normalized by total energy in original signal
#        energy_z_x = sum(np.abs(gz_x[:-p])**2)/sum(np.abs(x_eqdist)**2)*10**6
#        energy_z_y = sum(np.abs(gz_y[:-p])**2)/sum(np.abs(y_eqdist)**2)*10**6
#        Eg_z_x.append((ftype,energy_z_x))
#        Eg_z_y.append((ftype,energy_z_y))
#   
#        ## plot
#        #plt.close('all')
#        #np.set_printoptions(precision=2)
#        #fig = plt.figure()
#        #ax_lin,ax_lpe = fig.add_subplot(211),fig.add_subplot(212)
#        #fig.text(0.99, 0.96,fname[:len(fname)-4],fontsize=10,color='red',va='baseline',ha='right',multialignment='left')
#        #if ftype=='healthy':
#        #    fig.text(0.32, 0.955, 'HEALTHY ('+clock_type+')',fontsize=15,color='black',va='baseline',ha='right',multialignment='left')
#        #elif ftype=='impaired':
#        #    fig.text(0.32, 0.955, 'IMPAIRED ('+clock_type+')',fontsize=15,color='black',va='baseline',ha='right',multialignment='left')
#        ## predict x[n] from z[n] = x[n]+j*y[n]
#        #ax_lin.set_title('Linear Prediction of x[n] from q[n]=x[n]+j*y[n]',fontsize=20)
#        #ax_lin.plot(x_eqdist,color='blue',lw=2,label=r'$x_{eqdist}[n]$')
#        #ax_lin.plot(zhat_x,'k-.',lw=3,label=r'$\hat{x}[n]$')
#        #ax_lin.legend(loc='best')
#        #ax_lpe.plot(gz_x[:-p],color='red',lw=2,label='$g[n]$')
#        #ax_lpe.legend(loc='best')
#        #ax_lpe.text(10,1.5,r'$E_g/E_{total} \times 10^2 = $'+str(np.round_(energy_z_x,3)),fontsize=15)
#        #ax_lpe.text(10,2.3,'p = '+str(p),fontsize=15)
#        ##ax_lin.set_ylim(bottom=-70,top=70)
#        ##ax_lpe.set_ylim(bottom=-3,top=3)
#        #fig.savefig(path+'figs_raw/'+fname[:len(fname)-4]+'/complex_lpc_x_p'+str(p)+'_'+clock_type+'_'+fname[:len(fname)-4]+'.png')
#        ## predict y[n] from z[n] = x[n]+j*y[n]
#        #ax_lpe.clear()
#        #ax_lin.set_title('Linear Prediction of y[n] from q[n]=x[n]+j*y[n]',fontsize=20)
#        #ax_lin.plot(y_eqdist,color='blue',lw=2,label=r'$y_{eqdist}[n]$')
#        #ax_lin.plot(zhat_y,'k-.',lw=3,label=r'$\hat{y}[n]$')
#        #ax_lin.legend(loc='best')
#        #ax_lpe.plot(gz_y[:-p],color='red',lw=2,label='$g[n]$')
#        #ax_lpe.legend(loc='best')
#        #ax_lpe.text(10,2.3,'p = '+str(p),fontsize=15)
#        #ax_lpe.text(10,1.5,r'$E_g/E_{total} \times 10^2 = $'+str(np.round_(energy_z_y,3)),fontsize=15)
#        ##ax_lin.set_ylim(bottom=-70,top=70)
#        ##ax_lpe.set_ylim(bottom=-3,top=3)
#        #fig.savefig(path+'figs_raw/'+fname[:len(fname)-4]+'/complex_lpc_y_p'+str(p)+'_'+clock_type+'_'+fname[:len(fname)-4]+'.png')
#        ## the whole clock from z[n] = x[n]+j*y[n]
#        #fig2 = plt.figure()
#        #ax2 = fig2.add_subplot(111)
#        #ax2.clear()
#        #ax2.plot(x_eqdist,y_eqdist,color='blue',lw=2,label='Drawn Circle')
#        #ax2.plot(zhat_x,zhat_y,'k-.',label='Lin. Prediction')
#        ##ax2.set_xlim(-5.5,5.5)
#        ##ax2.set_ylim(-5.5,5.5)
#        #ax2.legend(loc='best',frameon=False)
#        #fig2.savefig('complex_lpc_wholeClock_p'+str(p)+'_'+clock_type+'_'+fname[:len(fname)-4]+'.png')
#
## compare total energy in linear prediction error between the drawings of healthy
## and impaired patients
#plt.close('all')
#fig = plt.figure()
#ax = fig.add_subplot(111)
#
#error_x_healthy = [abs(elt[1]) for elt in Eg_z_x if elt[0]=='healthy']
#error_x_impaired = [abs(elt[1]) for elt in Eg_z_x if elt[0]=='impaired']
#error_y_healthy = [abs(elt[1]) for elt in Eg_z_y if elt[0]=='healthy']
#error_y_impaired = [abs(elt[1]) for elt in Eg_z_y if elt[0]=='impaired']
#
#np.savetxt('complex_lpc_error_x_healthy_'+clock_type+'.txt',error_x_healthy)
#np.savetxt('complex_lpc_error_x_impaired_'+clock_type+'.txt',error_x_impaired)
#np.savetxt('complex_lpc_error_y_healthy_'+clock_type+'.txt',error_y_healthy)
#np.savetxt('complex_lpc_error_y_impaired_'+clock_type+'.txt',error_y_impaired)

make_hists(clock_type)
