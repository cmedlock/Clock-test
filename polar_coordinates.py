# -*- coding: utf-8 -*-
# plot the circles in polar coordinates
import math
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
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

    x,y,t = np.array(x)-x[0],np.array(y)-y[0],np.array(t)-t[0]
    
    # compensate for non-constant velocity
    N_orig = len(x)

    # calculate average distance between points
    dists = []
    for w in range(1,len(x)):
        dx,dy = x[w]-x[w-1],y[w]-y[w-1]
        dist = math.sqrt(dx**2+dy**2)
        dists.append(dist)
    dist_avg = mean(dists)
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
        n_segments = ceil(dist/dist_avg)*100
        for r in range(1,int(n_segments)):
            x_new = x[w]+r*dx/n_segments
            y_new = y[w]+r*dy/n_segments
            x_interp.append(x_new)
            y_interp.append(y_new)
    x_interp.append(x[-1])
    y_interp.append(y[-1])
    #print '\naverage distance between interpolated points is ',dist_avg_interp
    #print 'total distance is now ',sum(dists_interp)

    # start from the first point and find the ones that are 
    # approximately a distance dist_avg from each other
    x_eqdist,y_eqdist = [x_interp[0]],[y_interp[0]]
    idx = 0
    for k in range(len(x)-1):
        dist_total = 0
        for j in range(idx,len(x_interp)-1):
            dx,dy = x_interp[j+1]-x_interp[j],y_interp[j+1]-y_interp[j]
            dist_total += math.sqrt(dx**2+dy**2)
            #if abs(dist_total-dist_avg)<0.01:
            if abs(dist_total-dist_avg)<dist_avg/100.:
                idx = j+1
	        break
        x_eqdist.append(x_interp[idx])
        y_eqdist.append(y_interp[idx])
    x_eqdist,y_eqdist = np.array(x_eqdist)-x_eqdist[0],np.array(y_eqdist)-y_eqdist[0]
    # normalize to unit energy, so can compare between drawings
    x_eqdist,y_eqdist = x_eqdist/sum(x_eqdist**2),y_eqdist/sum(y_eqdist**2)

    # find COM
    x_com = np.mean(x_eqdist)
    y_com = np.mean(y_eqdist)

    # get r and theta
    r,theta = [],[]
    for w in range(len(x)):
        dx,dy = x_eqdist[w]-x_com,y_eqdist[w]-y_com
        dist = sqrt(dx**2+dy**2)
        angle = math.atan2(dy,dx)
        r.append(dist)
        theta.append(angle)

    # figure settings
    linecolor_xy = 'blue'
    linecolor_dft = 'red'
    x_axis_fontsize = 15
    y_axis_fontsize = 20

    # initialize figures
    fig_xt,fig_yt,fig_xy,fig_rtheta = plt.figure(),plt.figure(),plt.figure(),plt.figure()
    fig_xt.subplots_adjust(hspace=0.6,left=0.15)
    fig_yt.subplots_adjust(hspace=0.6,left=0.15)
    # declare subplots
    xt,yt,xy,rtheta = fig_xt.add_subplot(411),fig_yt.add_subplot(411),fig_xy.add_subplot(111),fig_rtheta.add_subplot(111)
    xt_eqdist,yt_eqdist = fig_xt.add_subplot(412),fig_yt.add_subplot(412)
    dftxk,dftyk = fig_xt.add_subplot(413),fig_yt.add_subplot(413)
    dftxk_zoom,dftyk_zoom = fig_xt.add_subplot(414),fig_yt.add_subplot(414)
    # draw circles (actual and ideal)
    xy.plot(x,y,lw=2,label='$x[n],y[n]$')
    xy.plot(x_true,y_true,'k-.',lw=3,label='$x_{true}[n],y_{true}[n]$')
    xy.legend(loc='best',frameon=False,fontsize=20)
    # draw circles in polar coordinates
    rtheta.plot(theta,r,label=r'$\theta[n],r[n]$')
    rtheta.legend(loc='best',frameon=False,fontsize=20)
    # draw x[n]
    xt.plot(x,lw=2)
    xt_eqdist.plot(x_eqdist,lw=2,label=r'$x_{eqdist}[n]$')
    xt_eqdist.plot(x_true,'k-.',lw=3,label=r'$x_{true}[n]$')
    xt_eqdist.legend(loc='best',frameon=False)
    dftxk.plot(omega_posfreq,np.abs(dftx_posfreq),linecolor_dft,lw=2)
    dftxk_zoom.plot(omega_posfreq_zoom,np.abs(dftx_posfreq_zoom),linecolor_dft,lw=2)
    # draw y[n]
    yt.plot(y,lw=2)
    yt_eqdist.plot(y_eqdist,lw=2,label=r'$y_{eqdist}[n]$')
    yt_eqdist.plot(y_true,'k-.',lw=3,label=r'$y_{true}[n]$')
    yt_eqdist.legend(loc='best',frameon=False)
    dftyk.plot(omega_posfreq,np.abs(dfty_posfreq),linecolor_dft,lw=2)
    dftyk_zoom.plot(omega_posfreq_zoom,np.abs(dfty_posfreq_zoom),linecolor_dft,lw=2)

    # equalize axis scales for circle drawing
    if max(x)-min(x)>max(y)-min(y):
        ax_range = max(x)-min(x)+20
        xy.set_xlim(min(x)-10,max(x)+10)
        xy.set_ylim(min(y)-10,min(y)+ax_range)
    else:
        ax_range = max(y)-min(y)+20
        xy.set_xlim(min(x)-10,min(x)+ax_range)
        xy.set_ylim(min(y)-10,max(y)+10)
    plt.axis('equal')

    # set axis limits
    rtheta.set_xlim(left=-pi-0.05,right=pi+0.05)
    rtheta.set_ylim(bottom=min(r)*1.2,top=max(r)*1.2)
    
    xt.set_xlim(right=len(x))
    xt_eqdist.set_xlim(right=len(x_eqdist))
    dftxk.set_ylim(bottom=min(np.abs(dftx_posfreq[1:]))/10,top=max(np.abs(dftx))*10)
    dftxk.set_yscale('log')
    dftxk_zoom.set_xlim(right=max(omega_posfreq_zoom))
    dftxk_zoom.set_ylim(bottom=min(np.abs(dftx_posfreq_zoom[1:]))/10,top=max(np.abs(dftx_zoom))*10)
    dftxk_zoom.set_yscale('log')

    yt.set_xlim(right=len(y))
    yt_eqdist.set_xlim(right=len(y_eqdist))
    dftyk.set_ylim(bottom=min(np.abs(dfty_posfreq[1:]))/10,top=max(np.abs(dfty))*10)
    dftyk.set_yscale('log')
    dftyk_zoom.set_xlim(right=max(omega_posfreq_zoom))
    dftyk_zoom.set_ylim(bottom=min(np.abs(dfty_posfreq_zoom[1:]))/10,top=max(np.abs(dfty_zoom))*10)
    dftyk_zoom.set_yscale('log')

    # set axis labels
    xy.set_xlabel(r'$x$',fontsize=x_axis_fontsize)
    xy.set_ylabel(r'$y$',fontsize=y_axis_fontsize)
    
    rtheta.set_xlabel('$\theta$',fontsize=x_axis_fontsize)
    rtheta.set_ylabel('$r$',fontsize=y_axis_fontsize)

    xt.set_xlabel(r'$n$',fontsize=x_axis_fontsize)
    xt.set_ylabel(r'$x[n]$',color=linecolor_xy,fontsize=y_axis_fontsize)
    xt_eqdist.set_xlabel(r'$n$',fontsize=x_axis_fontsize)
    xt_eqdist.set_ylabel('$x_{eqdist}[n]$,\n$x_{true}[n]$',fontsize=y_axis_fontsize-3)
    dftxk.set_xlabel(r'$\omega/\pi$',fontsize=x_axis_fontsize)
    dftxk_zoom.set_xlabel(r'$\omega/\pi$',fontsize=x_axis_fontsize)
    for a1 in xt.get_yticklabels():
        a1.set_color(linecolor_xy)
    for a2 in xt_eqdist.get_yticklabels():
        a2.set_color(linecolor_xy)

    yt.set_xlabel(r'$n$',fontsize=x_axis_fontsize)
    yt.set_ylabel(r'$y[n]$',color=linecolor_xy,fontsize=y_axis_fontsize)
    yt_eqdist.set_xlabel(r'$n$',fontsize=x_axis_fontsize)
    yt_eqdist.set_ylabel('$y_{eqdist}[n]$,\n$y_{true}[n]$',fontsize=y_axis_fontsize-3)
    dftyk.set_xlabel(r'$\omega/\pi$',fontsize=x_axis_fontsize)
    dftyk_zoom.set_xlabel(r'$\omega/\pi$',fontsize=x_axis_fontsize)
    for b1 in yt.get_yticklabels():
        b1.set_color(linecolor_xy)
    for b2 in yt_eqdist.get_yticklabels():
        b2.set_color(linecolor_xy)
        
    dftxk.set_ylabel(r'$|X_{eqdist}[k]|$',color=linecolor_dft,fontsize=y_axis_fontsize)
    dftyk.set_ylabel(r'$|Y_{eqdist}[k]|$',color=linecolor_dft,fontsize=y_axis_fontsize)
    for c1 in dftxk.get_yticklabels():
        c1.set_color(linecolor_dft)
    for c2 in dftyk.get_yticklabels():
        c2.set_color(linecolor_dft)

    dftxk_zoom.set_ylabel(r'$|X_{eqdist}[k]|$',color=linecolor_dft,fontsize=y_axis_fontsize)
    dftyk_zoom.set_ylabel(r'$|Y_{eqdist}[k]|$',color=linecolor_dft,fontsize=y_axis_fontsize)
    for d1 in dftxk_zoom.get_yticklabels():
        d1.set_color(linecolor_dft)
    for d2 in dftyk_zoom.get_yticklabels():
        d2.set_color(linecolor_dft)
    
    # add drawing type (healthy or impaired) and file name
    fig_xy.text(0.99, 0.96,fname[:len(fname)-4],fontsize=10,color='red',va='baseline',ha='right',multialignment='left')
    fig_xt.text(0.99, 0.96,fname[:len(fname)-4],fontsize=10,color='red',va='baseline',ha='right',multialignment='left')
    fig_yt.text(0.99, 0.96,fname[:len(fname)-4],fontsize=10,color='red',va='baseline',ha='right',multialignment='left')

    if 'YDU' in fname:
        fig_xy.text(0.25, 0.955, 'HEALTHY',fontsize=15,color='black',va='baseline',ha='right',multialignment='left')
        fig_xt.text(0.25, 0.955, 'HEALTHY',fontsize=15,color='black',va='baseline',ha='right',multialignment='left')
        fig_yt.text(0.25, 0.955, 'HEALTHY',fontsize=15,color='black',va='baseline',ha='right',multialignment='left')
    elif 'CIN' in fname:
        fig_xy.text(0.25, 0.955, 'IMPAIRED',fontsize=15,color='black',va='baseline',ha='right',multialignment='left')
        fig_xt.text(0.25, 0.955, 'IMPAIRED',fontsize=15,color='black',va='baseline',ha='right',multialignment='left')
        fig_yt.text(0.25, 0.955, 'IMPAIRED',fontsize=15,color='black',va='baseline',ha='right',multialignment='left')
    else:
        print 'not a valid filename'

    # save figures
    fig_xy.savefig(path+'figs_raw/'+fname[:len(fname)-4]+'/xy_'+clock_type+'_'+fname[:len(fname)-4]+'.png')
    fig_rtheta.savefig(path+'figs_raw/'+fname[:len(fname)-4]+'/xy_polar_'+clock_type+'_'+fname[:len(fname)-4]+'.png')
    fig_xt.savefig(path+'figs_raw/'+fname[:len(fname)-4]+'/x_true_'+clock_type+'_'+fname[:len(fname)-4]+'.png')
    fig_yt.savefig(path+'figs_raw/'+fname[:len(fname)-4]+'/y_true_'+clock_type+'_'+fname[:len(fname)-4]+'.png')

    plt.close('all')