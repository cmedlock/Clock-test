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
clock_type = 'COMMAND'

# save interesting quantities to compare between healthy and impaired patients
# format is [[('healthy',ak_1 val),('healthy',ak_1 val),('impaired',ak_1 val),...]
#            [ same thing for ak_2 ],...
#            [ same thing for ak_p]]
Eg_x,Eg_y = [],[] # energy in g[n] (i.e. linear prediction error)

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

    # subtract mean values so there is no DC term adding an extra pole
    x_eqdist = [elt-mean(x_eqdist) for elt in x_eqdist]
    y_eqdist = [elt-mean(y_eqdist) for elt in y_eqdist]
    
    x_eqdist,y_eqdist = np.array(x_eqdist),np.array(y_eqdist)
    
    # linearly predict z[n] = x[n] + j*y[n] by multiplying each previous value
    # of z[n] by e^(j*omega), where omega = arctan(y[n]/x[n])
    # i.e. zhat[0] = x[0] + j*y[0] = xhat[0] + j*yhat[0]
    #      zhat[1] = (xhat[0] + j*yhat[0]) * e^(j*arctan(yhat[0]/xhat[0]))
    #      etc.
    z = x_eqdist+1j*y_eqdist
    omega = 2*math.pi/float(len(x_eqdist))
    #z = np.exp(1j*2.*math.pi/251.*np.arange(251))
    inst_freq = np.angle(z)
    zhat = [z[0]]
    for w in range(1,len(z)):
        zhat.append(np.exp(1j*omega)*zhat[w-1])
    # linear predictions of x and y separately
    xhat,yhat = np.real(zhat),np.imag(zhat)
    ## get the phase shift right using matched filtering
    min_dist = 10**7
    best_shift = 0
    for shift in range(len(xhat)):
        xhat_shift = np.concatenate((xhat[shift:],xhat[:shift]))
        yhat_shift = np.concatenate((yhat[shift:],yhat[:shift]))
        dist = 0
        for w in range(len(xhat)):
            dx = xhat_shift[w]-x_eqdist[w]
            dy = yhat_shift[w]-y_eqdist[w]
            dist += (dx**2+dy**2)**0.5
        if dist<min_dist:
            min_dist = dist
            best_shift = shift
    xhat = np.concatenate((xhat[best_shift:],xhat[:best_shift]))
    yhat = np.concatenate((yhat[best_shift:],yhat[:best_shift]))
    # approximate amplitude by using average distance from origin
    dist_eqdist,dist_hat = 0,0
    for w in range(len(x_eqdist)):
        dist_eqdist += (x_eqdist[w]**2+y_eqdist[w]**2)**0.5
        dist_hat += (xhat[w]**2+yhat[w]**2)**0.5
    avg_dist_eqdist = dist_eqdist/float(len(x_eqdist))
    avg_dist_hat = dist_hat/float(len(xhat))
    xhat = xhat*avg_dist_eqdist/avg_dist_hat
    yhat = yhat*avg_dist_eqdist/avg_dist_hat
    
    # linear prediction error
    gx = xhat-x_eqdist
    gy = yhat-y_eqdist
    energy_x = sum(gx**2)/sum(x_eqdist**2)
    energy_y = sum(gy**2)/sum(y_eqdist**2)
    
    # store the coefficients for comparison between the drawings of healthy
    # and impaired patients
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
    ax_lin.set_title('Linear Prediction: x[n]',fontsize=20)
    ax_lin.plot(x_eqdist,color='blue',lw=2,label=r'$x_{eqdist}[n]$')
    ax_lin.plot(xhat,'k-.',lw=3,label=r'$\hat{x}[n]$')
    ax_lin.legend(loc='best')
    ax_lpe.plot(gx,color='red',lw=2,label='$g[n]$')
    ax_lpe.legend(loc='best')
    ax_lpe.text(0.05,0.15,r'$E_g/E_{total}= $'+str(np.round_(energy_x,3)),transform=ax_lpe.transAxes,fontsize=15,verticalalignment='top')
    fig.savefig(path+'figs_raw/'+fname[:len(fname)-4]+'/complexlpc_x_'+clock_type+'_'+fname[:len(fname)-4]+'.png')

    ax_lin.clear()
    ax_lpe.clear()
    # predict y[n]
    ax_lin.set_title('Linear Prediction: y[n]',fontsize=20)
    ax_lin.plot(y_eqdist,color='blue',lw=2,label=r'$y_{eqdist}[n]$')
    ax_lin.plot(yhat,'k-.',lw=3,label=r'$\hat{y}[n]$')
    ax_lin.legend(loc='best')
    ax_lpe.plot(gx,color='red',lw=2,label='$g[n]$')
    ax_lpe.legend(loc='best')
    ax_lpe.text(0.05,0.15,r'$E_g/E_{total}= $'+str(np.round_(energy_y,3)),transform=ax_lpe.transAxes,fontsize=15,verticalalignment='top')
    fig.savefig(path+'figs_raw/'+fname[:len(fname)-4]+'/complexlpc_y_'+clock_type+'_'+fname[:len(fname)-4]+'.png')
    
    # the whole clock
    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111)
    ax2.plot(x_eqdist,y_eqdist,color='blue',lw=2,label='True Clock')
    ax2.plot(xhat,yhat,'k-.',label='Lin. Prediction')
    ax2.legend(loc='best',frameon=False)
    fig2.savefig(path+'figs_raw/'+fname[:len(fname)-4]+'/complexlpc_wholeClock_'+clock_type+'_'+fname[:len(fname)-4]+'.png')
    
    # the instataneous frequency
    fig3 = plt.figure()
    ax3 = fig3.add_subplot(111)
    ax3.plot(inst_freq,color='blue',lw=2,label='Inst. Freq.')
    ax3.legend(loc='best',frameon=False)
    fig3.savefig(path+'figs_raw/'+fname[:len(fname)-4]+'/instfreq_'+clock_type+'_'+fname[:len(fname)-4]+'.png')

# compare total energy in linear prediction error 
# between the drawings of healthy and impaired patients
plt.close('all')

# energy in linear prediction error compared to energy of signal
ct.make_hist([elt[1] for elt in Eg_x if elt[0]=='healthy'],
             [elt[1] for elt in Eg_x if elt[0]=='impaired'],
             10,'Relative Energy in Linear Prediciton Error (Complex)','complexEg_x_'+clock_type,path)
ct.make_hist([elt[1] for elt in Eg_y if elt[0]=='healthy'],
             [elt[1] for elt in Eg_y if elt[0]=='impaired'],
             10,'Relative Energy in Linear Prediciton Error (Complex)','complexEg_y_'+clock_type,path)
