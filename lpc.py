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

# save interesting quantities
# model order
p = 2
ak_x_coeffs,ak_y_coeffs = [],[]
ck_x_coeffs,ck_y_coeffs = [],[]
for w in range(p):
    ak_x_coeffs.append([])
    ak_y_coeffs.append([])
    ck_x_coeffs.append([])
    ck_y_coeffs.append([])
    
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
    if not os.path.exists(path+'figs_raw/'+fname[:len(fname)-4]+'/LPC'):
        os.makedirs(path+'figs_raw/'+fname[:len(fname)-4]+'/LPC')

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

    # form all-pole model using Y-W eqns

    p = 2 # model order
    # 'circular' autocorrelation
    rxx,ryy = [],[]
    x_periodic,y_periodic = np.concatenate((x_eqdist,x_eqdist)),np.concatenate((y_eqdist,y_eqdist))
    for w in range(p+1):
        rxx.append(np.dot(x_eqdist,x_periodic[w:w+len(x_eqdist)]))
        ryy.append(np.dot(y_eqdist,y_periodic[w:w+len(y_eqdist)]))
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
    # LPC cepstrum
    ck_x,ck_y = [ak_x[0]],[ak_y[0]]
    for k in range(2,p+1):
        x1,y1 = ak_x[k-1],ak_y[k-1]
        for m in range(1,k):
            x1 += float(m)/float(k)*ak_x[m-1]*ck_x[k-m-1]
            y1 += float(m)/float(k)*ak_y[m-1]*ck_y[k-m-1]
        ck_x.append(x1)
        ck_y.append(y1)

    # store the coefficients for comparison between the drawings of healthy
    # and impaired patients
    for m in range(p):
        ak_x_coeffs[m].append((ftype,ak_x[m]))
        ak_y_coeffs[m].append((ftype,ak_y[m]))
        ck_x_coeffs[m].append((ftype,ck_x[m]))
        ck_y_coeffs[m].append((ftype,ck_y[m]))
        
    # plot
    plt.close('all')
    fig_x = plt.figure()
    ax_xspec,ax_xcep = fig_x.add_subplot(211),fig_x.add_subplot(212)
    ax_xspec.stem(ak_x)
    ax_xcep.stem(ck_x)
    ax_xspec.set_xlim(left=-1,right=len(ak_x))
    ax_xcep.set_xlim(left=-1,right=len(ck_x))
    ymin_spec = min(ak_x)*1.2 if min(ak_x)<0 else 0
    ymin_cep = min(ck_x)*1.2 if min(ck_x)<0 else 0
    ax_xspec.set_ylim(top=max(ak_x)*1.2,bottom=ymin_spec)
    ax_xcep.set_ylim(top=max(ck_x)*1.2,bottom=ymin_cep)
    ax_xspec.set_ylabel('$a_{kx}$',fontsize=20)
    ax_xcep.set_ylabel('$c_{kx}$',fontsize=20)
    
    fig_y = plt.figure()
    ax_yspec,ax_ycep = fig_y.add_subplot(211),fig_y.add_subplot(212)
    ax_yspec.stem(ak_y)
    ax_ycep.stem(ck_y)
    ax_yspec.set_xlim(left=-1,right=len(ak_y))
    ax_ycep.set_xlim(left=-1,right=len(ck_y))
    ymin_spec = min(ak_y)*1.2 if min(ak_y)<0 else 0
    ymin_cep = min(ck_y)*1.2 if min(ck_y)<0 else 0
    ax_yspec.set_ylim(top=max(ak_y)*1.2,bottom=ymin_spec)
    ax_ycep.set_ylim(top=max(ck_y)*1.2,bottom=ymin_cep)
    ax_yspec.set_ylabel('$a_{ky}$',fontsize=20)
    ax_ycep.set_ylabel('$c_{ky}$',fontsize=20)

    fig_x.savefig(path+'figs_raw/'+fname[:len(fname)-4]+'/LPC/lpc_x_'+clock_type+'_'+fname[:len(fname)-4]+'.png')
    fig_y.savefig(path+'figs_raw/'+fname[:len(fname)-4]+'/LPC/lpc_y_'+clock_type+'_'+fname[:len(fname)-4]+'.png')

# compare LPC spectrum and cepstrum for the drawings of healthy vs. impaired patients
for w in range(p):
    # LPC spectrum
    ct.make_hist([elt[1] for elt in ak_x_coeffs[w] if elt[0]=='healthy'],
                [elt[1] for elt in ak_x_coeffs[w] if elt[0]=='impaired'],
                10,clock_type+' Clock LPC: a_'+str(w)+' (x)','a'+str(w)+'_x_'+clock_type,path)
    ct.make_hist([elt[1] for elt in ak_y_coeffs[w] if elt[0]=='healthy'],
                [elt[1] for elt in ak_y_coeffs[w] if elt[0]=='impaired'],
                10,clock_type+' Clock LPC: a_'+str(w)+' (y)','a'+str(w)+'_y_'+clock_type,path)
    # LPC cepstrum
    ct.make_hist([elt[1] for elt in ck_x_coeffs[w] if elt[0]=='healthy'],
                [elt[1] for elt in ck_x_coeffs[w] if elt[0]=='impaired'],
                10,clock_type+' Clock LPC Cepstrum: c_'+str(w)+' (x)','c'+str(w)+'_x_'+clock_type,path)
    ct.make_hist([elt[1] for elt in ak_y_coeffs[w] if elt[0]=='healthy'],
                [elt[1] for elt in ak_y_coeffs[w] if elt[0]=='impaired'],
                10,clock_type+' Clock LPC Cepstrum: c_'+str(w)+' (y)','c'+str(w)+'_y_'+clock_type,path)
