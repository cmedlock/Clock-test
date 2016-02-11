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

# save interesting quantities
corr_x,corr_y = [],[]
Ediff_x,Ediff_y = [],[]

for fname in dirs[:2]:
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

    # copy or command clock?
    clock_type = 'COPY'
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

    # step 1: normalize the size of the word such that all coordinates
    # are linearly positioned between 0 and 127
    x = 127*(x-min(x))/(max(x)-min(x))
    y = 127*(y-min(y))/(max(y)-min(y))
    
    # step 2: resample to compensate for different velocities and divide
    # the symbol into frames of length L = 50 samples offset by R = 30 samples
    # from each other
    N_orig = len(x)
    N_new = 250

    # calculate total and average distance between points
    dists = []
    for w in range(1,len(x)):
        dx,dy = x[w]-x[w-1],y[w]-y[w-1]
        dist = math.sqrt(dx**2+dy**2)
        dists.append(dist)
    dist_avg = mean(dists)
    dist_total = sum(dists)
    
    #### now want to get N_orig evenly-spaced points along the curve
    # now want to get 250 evenly-spaced points along the curve

    # generate a much longer array with 199 linearly-interpolated 
    # points between the actual data points
    x_interp,y_interp = [],[]
    for w in range(len(x)-1):
        x_interp.append(x[w])
        y_interp.append(y[w])
        dx,dy = x[w+1]-x[w],y[w+1]-y[w]
        dist = math.sqrt(dx**2+dy**2)
        n_segments = ceil(dist/dist_avg)*100
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
    #for k in range(len(x)-1):
    for k in range(N_new):
        dist_sofar = 0
        for j in range(idx,len(x_interp)-1):
            dx,dy = x_interp[j+1]-x_interp[j],y_interp[j+1]-y_interp[j]
            dist_sofar += math.sqrt(dx**2+dy**2)
            if abs(dist_sofar-dist_total/250.)<dist_total/(250.*100.):
                idx = j+1
	        break
        x_eqdist.append(x_interp[idx])
        y_eqdist.append(y_interp[idx])
    x_eqdist,y_eqdist = np.array(x_eqdist)-x_eqdist[0],np.array(y_eqdist)-y_eqdist[0]
    
    # get frames
    L = 50 # frame size
    R = 30 # frame shift
    F = (N_new-L)/R+1 # total number of frames
    
    xframes,yframes = [],[]
    for w in range(F):
        xframes.append(x_eqdist[w*R:w*R+L])
        yframes.append(y_eqdist[w*R:w*R+L])
    
    # step 3: window each frame to minimize spectral leakage
    n = np.arange(L)
    window = 0.54-0.46*np.cos(2*pi*n/(L-1))
    xframes = [np.array(elt)*window for elt in xframes]
    yframes = [np.array(elt)*window for elt in yframes]
    
    # step 4: calculate the LPC cepstrum of each frame
    for w in range(F):
        xframe,yframe = xframes[w],yframes[w]
        # calculate deterministic autocorrelation
        rxx,ryy = [],[]
        for r in range(2*len(xframe)-1):
            xval,yval = 0,0
            for d in range(len(xframe)):
   	        if -1<len(xframe)-1-r+d<len(xframe):
                    xval += xframe[d]*xframe[len(xframe)-1-r+d]
                    yval += yframe[d]*yframe[len(yframe)-1-r+d]
            rxx.append(xval)
            ryy.append(yval)
        # model order
        p = 12
        center_idx = len(xframe)-1
        D_x,D_y = np.array(rxx[center_idx+1:center_idx+1+p]),np.array(ryy[center_idx+1:center_idx+1+p])
        W_x,W_y = np.empty((p,p)),np.empty((p,p))
        for row in range(p):
            for column in range(row,p):
                W_x[row][column] = rxx[center_idx+column-row]
                W_x[column][row] = rxx[center_idx+column-row]
                W_y[row][column] = ryy[center_idx+column-row]
                W_y[column][row] = ryy[center_idx+column-row]
        # LPC spectrum
        W_x_inv,W_y_inv = np.linalg.inv(W_x),np.linalg.inv(W_y)
        ak_x,ak_y = np.dot(W_x_inv,D_x),np.dot(W_y_inv,D_y)
        # LPC cepstrum
        ck_x,ck_y = [ak_x[0]],[ak_y[0]]
        for k in range(2,p+1):
            x1,y1 = ak_x[k-1],ak_y[k-1]
            for m in range(1,k):
                x1 += m/k*ak_x[m-1]*ck_x[k-m-1]
                y1 += m/k*ak_y[m-1]*ck_y[k-m-1]
            ck_x.append(x1)
            ck_y.append(y1)
        # plot
        plt.close('all')
        fig = plt.figure()
        ax_x,ax_y = fig.add_subplot(211),fig.add_subplot(212)
        ax_x.stem(ck_x)
        ax_y.stem(ck_y)
        ax_x.set_xlim(left=-1)
        ax_y.set_xlim(left=-1)
        ax_x.set_ylim(top=max(ck_x)*1.2)
        ax_y.set_ylim(top=max(ck_y)*1.2)
        ax_x.set_ylabel('$c_{kx}$')
        ax_y.set_ylabel('$c_{ky}$')
        fig1 = plt.figure()
        ax_x1,ax_y1 = fig1.add_subplot(211),fig1.add_subplot(212)
        ax_x1.stem(ak_x)
        ax_y1.stem(ak_y)
        ax_x1.set_xlim(left=-1)
        ax_y1.set_xlim(left=-1)
        ax_x1.set_ylim(top=max(ak_x)*1.2)
        ax_y1.set_ylim(top=max(ak_y)*1.2)
        ax_x1.set_ylabel('$a_{kx}$')
        ax_y1.set_ylabel('$a_{ky}$')

        fig.savefig(path+'figs_raw/'+fname[:len(fname)-4]+'/lpc_cepstrum_'+clock_type+'_'+fname[:len(fname)-4]+'.png')
        fig1.savefig(path+'figs_raw/'+fname[:len(fname)-4]+'/lpc_spectrum_'+clock_type+'_'+fname[:len(fname)-4]+'.png')
