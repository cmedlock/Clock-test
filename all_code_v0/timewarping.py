# -*- coding: utf-8 -*-
# 1.normalize for non-constant velocity
# 2.estimate ideal underlying sinusoids
# 3.compute correlation and mean squared difference
# between ideal sinusoids and actual drawing
# 4.draw clockface in polar coordinates
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

# path to output
if not os.path.exists(path+'figs_raw'):
    os.makedirs(path+'figs_raw')

# arrays to store the polar coordinates of all clocks
# each clock should be labeled as either healthy or impaired
x_all,y_all = [],[]
r_all,theta_all = [],[]

# save interesting quantities
corr_x,corr_y = [],[]
Ediff_x,Ediff_y = [],[]

# copy or command clock?
clock_type = 'COPY'

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
    #for k in range(len(x)-1):
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
    
    # store, plot outside of loop
    x_all.append((ftype,x_eqdist))
    y_all.append((ftype,y_eqdist))

#    # now want to estimate the frequency of the underlying sinusoid
#    # subtract mean values (zm = zero mean)
#    x_eqdist_zm,y_eqdist_zm = x_eqdist-mean(x_eqdist),y_eqdist-mean(y_eqdist)
#    
#    # DFS coefficients
#    dft_size = len(x_eqdist_zm)
#    k = np.arange(dft_size)
#    omega = 2.*k/float(dft_size)
#    dftx,dfty = np.fft.fft(x_eqdist_zm,n=dft_size),np.fft.fft(y_eqdist_zm,n=dft_size)
#
#    # k_near_pi is the k value for which w_k = 2*pi*k/N is closest to,
#    # but not larger than, pi
#    k_near_pi = 0
#    if dft_size%2==0:
#        k_near_pi = dft_size/2+1
#    else:
#        k_near_pi = math.ceil(dft_size/2)
#
#    # only use the positive frequencies
#    posfreq = k[:k_near_pi]
#    omega_posfreq = omega[:k_near_pi]
#    dftx_posfreq,dfty_posfreq = dftx[:k_near_pi],dfty[:k_near_pi]
#    #freq = 2*np.pi/dft_size*k
#    #freq = np.fft.fftfreq(n=1024,d=1/(2*np.pi)) # NB: centered around 0
#    
#    # take the frequency of the largest DFS coefficient to be
#    # the approximate frequency of the underlying sinusoid
#    abs_dftx_posfreq,abs_dfty_posfreq = list(np.abs(dftx_posfreq)),list(np.abs(dfty_posfreq))
#    k_true_x = abs_dftx_posfreq.index(max(abs_dftx_posfreq))
#    k_true_y = abs_dfty_posfreq.index(max(abs_dfty_posfreq))
#    w_true_x = 2*pi*k_true_x/dft_size
#    w_true_y = 2*pi*k_true_y/dft_size
#    #print '\nw_true_x = ',w_true_x,' and w_true_y = ',w_true_y
#    #print '\n',np.angle(dftx[1]),np.angle(dftx[-1]),' ---> ',np.angle(dftx[1])-np.angle(dftx[-1])
#    #print np.angle(dfty[1]),np.angle(dfty[-1]),' ---> ',np.angle(dfty[1])-np.angle(dfty[-1])
#
#    # calculate ideal underlying sinusoid
#    n = np.arange(len(x_eqdist))
#    x_true,y_true = np.cos(w_true_x*n),np.sin(w_true_y*n)
#    # use maximum correlation to determine phase and range of x or y to determine amplitude
#    phase_x,max_corr_x = 0,0
#    phase_y,max_corr_y = 0,0
#    for w in range(len(x_true)):
#        x_shifted = np.concatenate((x_true[w:],x_true[:w]))
#        y_shifted = np.concatenate((y_true[w:],y_true[:w]))
#        rho_x = np.dot(x_eqdist,x_shifted)
#        rho_y = np.dot(y_eqdist,y_shifted)
#        if rho_x>max_corr_x:
#            max_corr_x = rho_x
#            phase_x = w
#        if rho_y>max_corr_y:
#            max_corr_y = rho_y
#            phase_y = w
#    x_true = np.concatenate((x_true[phase_x:],x_true[:phase_x]))
#    y_true = np.concatenate((y_true[phase_y:],y_true[:phase_y]))
#    # use range of x or y to estimate amplitude
#    amp_x = (max(x_eqdist)-min(x_eqdist))/2
#    amp_y = (max(y_eqdist)-min(y_eqdist))/2
#    x_true = amp_x*x_true
#    y_true = amp_y*y_true
#    x_true,y_true = x_true+mean(x_eqdist),y_true+mean(y_eqdist)
#
#    # compute and save correlation between the two
#    corr_x.append(np.dot(x_eqdist,x_true))
#    corr_y.append(np.dot(y_eqdist,y_true))
#    
#    # compute and save mean squared difference between the two
#    Etrue_x,Etrue_y = sum(x_true**2),sum(y_true**2)
#    Ediff_x.append(sum((x_eqdist-x_true)**2))
#    Ediff_y.append(sum((y_eqdist-y_true)**2))
#
#    # zoom in on the DFS coefficients between ±(pi/8)
#
#    # extract the desired frequency band from x_eqdist[n],
#    # then recalculate the DFS coefficients
#        
#    # convolve (truncated) sinc with x_eqdist[n]
#    omega_c = pi/8
#    length_of_sinc = 601
#    x1,y1 = [],[]
#    for n in range(len(x_eqdist_zm)+length_of_sinc-1):
#        x1_of_n,y1_of_n = 0,0
#        for m in range(len(x_eqdist_zm)):
#            x1_of_n += x_eqdist_zm[m]*sinc(omega_c,n-(length_of_sinc-1)/2-m,length_of_sinc)
#            y1_of_n += y_eqdist_zm[m]*sinc(omega_c,n-(length_of_sinc-1)/2-m,length_of_sinc)
#        x1.append(x1_of_n)
#        y1.append(y1_of_n)
#    # check dft of x1
#    dftx1 = np.fft.fft(x1,n=len(x1))
#    #print 'check for x1: ',sum(x1),' = ',np.abs(dftx1[0]),' = ',8*np.abs(dftx[0]),'?'
#    
#    # decimate by 8 (d = decimated)
#    x1_d,y1_d = [],[]
#    for w in range(len(x1)):
#        if w%8==0:
#            x1_d.append(x1[w])
#            y1_d.append(y1[w])
#    
#    x1_d,y1_d = np.array(x1_d),np.array(y1_d)
#    
#    # redo the N-point DFT
#    dft_zoom_size = len(x1_d)
#    k_zoom = np.arange(dft_zoom_size)
#    omega_zoom = 2.*k_zoom/(8*float(dft_zoom_size)) # NB: these frequencies correspond to samples of X_eqdist(e^jw), not of X1D(e^jw)
#    dftx_zoom,dfty_zoom = np.fft.fft(x1_d,n=dft_zoom_size),np.fft.fft(y1_d,n=dft_zoom_size)
#    #print 'check 1: ',sum(x1_d),' = ',np.abs(dftx_zoom[0]),' = ',np.abs(dftx[0]),'?'
#    #print 'check 2: ',sum(y1_d),' = ',np.abs(dfty_zoom[0]),' = ',np.abs(dfty[0]),'?'
#
#    # k_zoom_near_pi is the smallest k value for which w_k = 2*pi*k/N is
#    # greater than pi
#    k_zoom_near_pi = 0
#    if dft_zoom_size%2==0:
#        k_zoom_near_pi = dft_zoom_size/2+1
#    else:
#        k_zoom_near_pi = math.ceil(dft_zoom_size/2)
#
#    # only use the positive frequencies for plotting
#    posfreq_zoom = k_zoom[:k_zoom_near_pi]
#    omega_posfreq_zoom = omega_zoom[:k_zoom_near_pi]
#    dftx_posfreq_zoom,dfty_posfreq_zoom = dftx_zoom[:k_zoom_near_pi],dfty_zoom[:k_zoom_near_pi]

    # plot circle in polar coordinates
    # find COM
    x_com = np.mean(x_eqdist)
    y_com = np.mean(y_eqdist)

    # get r and theta
    r,theta = [],[]
    for w in range(len(x_eqdist)):
        dx,dy = x_eqdist[w]-x_com,y_eqdist[w]-y_com
        dist = sqrt(dx**2+dy**2)
        angle = math.atan2(dy,dx)
        if angle<0:
            angle = angle+2*pi
        r.append(dist)
        theta.append(angle)
    r,theta = np.array(r),np.array(theta)
    # store, plot outside the loop
    r_all.append((ftype,r))
    theta_all.append((ftype,theta))

#    # figure settings
#    linecolor_xy = 'blue'
#    linecolor_dft = 'red'
#    x_axis_fontsize = 15
#    y_axis_fontsize = 20
#
#    # initialize figures
#    fig_xt,fig_yt,fig_xy,fig_rtheta = plt.figure(),plt.figure(),plt.figure(),plt.figure()
#    fig_xt.subplots_adjust(hspace=0.6,left=0.15)
#    fig_yt.subplots_adjust(hspace=0.6,left=0.15)
#    # declare subplots
#    xt,yt,xy,rtheta = fig_xt.add_subplot(411),fig_yt.add_subplot(411),fig_xy.add_subplot(111,aspect=1.0),fig_rtheta.add_subplot(111)
#    xt_eqdist,yt_eqdist = fig_xt.add_subplot(412),fig_yt.add_subplot(412)
#    dftxk,dftyk = fig_xt.add_subplot(413),fig_yt.add_subplot(413)
#    dftxk_zoom,dftyk_zoom = fig_xt.add_subplot(414),fig_yt.add_subplot(414)
#    # draw circles (actual and ideal)
#    xy.plot(x,y,lw=2,label='$x[n],y[n]$')
#    xy.plot(x_true,y_true,'k-.',lw=3,label='$x_{true}[n],y_{true}[n]$')
#    xy.legend(loc='best',frameon=False,fontsize=20)
#    # draw circles in polar coordinates
#    rtheta.plot(theta,r,label=r'$\theta[n],r[n]$')
#    rtheta.legend(loc='best',frameon=False,fontsize=20)
#    # draw x[n]
#    xt.plot(x,lw=2)
#    xt_eqdist.plot(x_eqdist,lw=2,label=r'$x_{eqdist}[n]$')
#    xt_eqdist.plot(x_true,'k-.',lw=3,label=r'$x_{true}[n]$')
#    xt_eqdist.legend(loc='best',frameon=False)
#    dftxk.plot(omega_posfreq,np.abs(dftx_posfreq),linecolor_dft,lw=2)
#    dftxk_zoom.plot(omega_posfreq_zoom,np.abs(dftx_posfreq_zoom),linecolor_dft,lw=2)
#    # draw y[n]
#    yt.plot(y,lw=2)
#    yt_eqdist.plot(y_eqdist,lw=2,label=r'$y_{eqdist}[n]$')
#    yt_eqdist.plot(y_true,'k-.',lw=3,label=r'$y_{true}[n]$')
#    yt_eqdist.legend(loc='best',frameon=False)
#    dftyk.plot(omega_posfreq,np.abs(dfty_posfreq),linecolor_dft,lw=2)
#    dftyk_zoom.plot(omega_posfreq_zoom,np.abs(dfty_posfreq_zoom),linecolor_dft,lw=2)
#
#    # set axis limits
#    xy.set_xlim(left=-10,right=140)
#    xy.set_ylim(bottom=-10,top=140)
#    
#    rtheta.set_xlim(left=-3.2,right=3.2)
#    rtheta.set_ylim(bottom=45,top=75)
#
#    xt.set_xlim(right=len(x))
#    xt_eqdist.set_xlim(right=len(x_eqdist))
#    dftxk.set_ylim(bottom=min(np.abs(dftx_posfreq[1:]))/10,top=max(np.abs(dftx))*10)
#    dftxk.set_yscale('log')
#    dftxk_zoom.set_xlim(right=max(omega_posfreq_zoom))
#    dftxk_zoom.set_ylim(bottom=min(np.abs(dftx_posfreq_zoom[1:]))/10,top=max(np.abs(dftx_zoom))*10)
#    dftxk_zoom.set_yscale('log')
#
#    yt.set_xlim(right=len(y))
#    yt_eqdist.set_xlim(right=len(y_eqdist))
#    dftyk.set_ylim(bottom=min(np.abs(dfty_posfreq[1:]))/10,top=max(np.abs(dfty))*10)
#    dftyk.set_yscale('log')
#    dftyk_zoom.set_xlim(right=max(omega_posfreq_zoom))
#    dftyk_zoom.set_ylim(bottom=min(np.abs(dfty_posfreq_zoom[1:]))/10,top=max(np.abs(dfty_zoom))*10)
#    dftyk_zoom.set_yscale('log')
#
#    # set axis labels
#    xy.set_xlabel(r'$y$',fontsize=x_axis_fontsize)
#    xy.set_ylabel(r'$x$',fontsize=y_axis_fontsize)
#
#    rtheta.set_xlabel(r'$\theta$',fontsize=x_axis_fontsize)
#    rtheta.set_ylabel(r'$r$',fontsize=y_axis_fontsize)
#
#    xt.set_xlabel(r'$n$',fontsize=x_axis_fontsize)
#    xt.set_ylabel(r'$x[n]$',color=linecolor_xy,fontsize=y_axis_fontsize)
#    xt_eqdist.set_xlabel(r'$n$',fontsize=x_axis_fontsize)
#    xt_eqdist.set_ylabel('$x_{eqdist}[n]$,\n$x_{true}[n]$',fontsize=y_axis_fontsize-3)
#    dftxk.set_xlabel(r'$\omega/\pi$',fontsize=x_axis_fontsize)
#    dftxk_zoom.set_xlabel(r'$\omega/\pi$',fontsize=x_axis_fontsize)
#    for a1 in xt.get_yticklabels():
#        a1.set_color(linecolor_xy)
#    for a2 in xt_eqdist.get_yticklabels():
#        a2.set_color(linecolor_xy)
#
#    yt.set_xlabel(r'$n$',fontsize=x_axis_fontsize)
#    yt.set_ylabel(r'$y[n]$',color=linecolor_xy,fontsize=y_axis_fontsize)
#    yt_eqdist.set_xlabel(r'$n$',fontsize=x_axis_fontsize)
#    yt_eqdist.set_ylabel('$y_{eqdist}[n]$,\n$y_{true}[n]$',fontsize=y_axis_fontsize-3)
#    dftyk.set_xlabel(r'$\omega/\pi$',fontsize=x_axis_fontsize)
#    dftyk_zoom.set_xlabel(r'$\omega/\pi$',fontsize=x_axis_fontsize)
#    for b1 in yt.get_yticklabels():
#        b1.set_color(linecolor_xy)
#    for b2 in yt_eqdist.get_yticklabels():
#        b2.set_color(linecolor_xy)
#        
#    dftxk.set_ylabel(r'$|X_{eqdist}[k]|$',color=linecolor_dft,fontsize=y_axis_fontsize)
#    dftyk.set_ylabel(r'$|Y_{eqdist}[k]|$',color=linecolor_dft,fontsize=y_axis_fontsize)
#    for c1 in dftxk.get_yticklabels():
#        c1.set_color(linecolor_dft)
#    for c2 in dftyk.get_yticklabels():
#        c2.set_color(linecolor_dft)
#
#    dftxk_zoom.set_ylabel(r'$|X_{eqdist}[k]|$',color=linecolor_dft,fontsize=y_axis_fontsize)
#    dftyk_zoom.set_ylabel(r'$|Y_{eqdist}[k]|$',color=linecolor_dft,fontsize=y_axis_fontsize)
#    for d1 in dftxk_zoom.get_yticklabels():
#        d1.set_color(linecolor_dft)
#    for d2 in dftyk_zoom.get_yticklabels():
#        d2.set_color(linecolor_dft)
#    
#    # add drawing type (healthy or impaired) and file name
#    fig_xy.text(0.99, 0.96,fname[:len(fname)-4],fontsize=10,color='red',va='baseline',ha='right',multialignment='left')
#    fig_rtheta.text(0.99, 0.96,fname[:len(fname)-4],fontsize=10,color='red',va='baseline',ha='right',multialignment='left')
#    fig_xt.text(0.99, 0.96,fname[:len(fname)-4],fontsize=10,color='red',va='baseline',ha='right',multialignment='left')
#    fig_yt.text(0.99, 0.96,fname[:len(fname)-4],fontsize=10,color='red',va='baseline',ha='right',multialignment='left')
#
#    if 'YDU' in fname:
#        fig_xy.text(0.32, 0.955, 'HEALTHY ('+clock_type+')',fontsize=15,color='black',va='baseline',ha='right',multialignment='left')
#        fig_rtheta.text(0.32, 0.955, 'HEALTHY ('+clock_type+')',fontsize=15,color='black',va='baseline',ha='right',multialignment='left')
#        fig_xt.text(0.32, 0.955, 'HEALTHY ('+clock_type+')',fontsize=15,color='black',va='baseline',ha='right',multialignment='left')
#        fig_yt.text(0.32, 0.955, 'HEALTHY ('+clock_type+')',fontsize=15,color='black',va='baseline',ha='right',multialignment='left')
#    elif 'CIN' in fname:
#        fig_xy.text(0.32, 0.955, 'IMPAIRED ('+clock_type+')',fontsize=15,color='black',va='baseline',ha='right',multialignment='left')
#        fig_rtheta.text(0.32, 0.955, 'IMPAIRED ('+clock_type+')',fontsize=15,color='black',va='baseline',ha='right',multialignment='left')
#        fig_xt.text(0.32, 0.955, 'IMPAIRED ('+clock_type+')',fontsize=15,color='black',va='baseline',ha='right',multialignment='left')
#        fig_yt.text(0.32, 0.955, 'IMPAIRED ('+clock_type+')',fontsize=15,color='black',va='baseline',ha='right',multialignment='left')
#    else:
#        print 'not a valid filename'
#
#    # save figures
#    fig_xy.savefig(path+'figs_raw/'+fname[:len(fname)-4]+'/xy_'+clock_type+'_'+fname[:len(fname)-4]+'.png')
#    fig_rtheta.savefig(path+'figs_raw/'+fname[:len(fname)-4]+'/xy_polar_'+clock_type+'_'+fname[:len(fname)-4]+'.png')
#    fig_xt.savefig(path+'figs_raw/'+fname[:len(fname)-4]+'/x_true_'+clock_type+'_'+fname[:len(fname)-4]+'.png')
#    fig_yt.savefig(path+'figs_raw/'+fname[:len(fname)-4]+'/y_true_'+clock_type+'_'+fname[:len(fname)-4]+'.png')

# separate polar coordinates into healthy and impaired
x_all_healthy = [elt[1] for elt in x_all if elt[0]=='healthy']
y_all_healthy = [elt[1] for elt in y_all if elt[0]=='healthy']
r_all_healthy = [elt[1] for elt in r_all if elt[0]=='healthy']
theta_all_healthy = [elt[1] for elt in theta_all if elt[0]=='healthy']
x_all_impaired = [elt[1] for elt in x_all if elt[0]=='impaired']
y_all_impaired = [elt[1] for elt in y_all if elt[0]=='impaired']
r_all_impaired = [elt[1] for elt in r_all if elt[0]=='impaired']
theta_all_impaired = [elt[1] for elt in theta_all if elt[0]=='impaired']

if not os.path.exists(path+'compare_healthy_impaired/all_clocks_overlaid'):
    os.makedirs(path+'compare_healthy_impaired/all_clocks_overlaid')

# plot
plt.close('all')
fig_healthy,fig_impaired = plt.figure(),plt.figure()
healthy,impaired = fig_healthy.add_subplot(111),fig_impaired.add_subplot(111)
fig_healthy.text(0.32, 0.955, 'HEALTHY ('+clock_type+')',fontsize=15,color='black',va='baseline',ha='right',multialignment='left')
fig_impaired.text(0.32, 0.955, 'IMPAIRED ('+clock_type+')',fontsize=15,color='black',va='baseline',ha='right',multialignment='left')

# just x
for w in range(len(x_all_healthy)):
    healthy.plot(x_all_healthy[w],color='blue')
for w in range(len(x_all_impaired)):
    impaired.plot(x_all_impaired[w],color='blue')
healthy.set_xlabel('n',fontsize=20)
healthy.set_ylabel('x',fontsize=20)
healthy.set_ylim(bottom=-20,top=140)
impaired.set_xlabel('n',fontsize=20)
impaired.set_ylabel('x',fontsize=20)
impaired.set_ylim(bottom=-20,top=140)

fig_healthy.savefig(path+'compare_healthy_impaired/all_clocks_overlaid/'+'/healthy_'+clock_type+'_clocks_x.png')
fig_impaired.savefig(path+'compare_healthy_impaired/all_clocks_overlaid/'+'/impaired_'+clock_type+'_clocks_x.png')

# just y
healthy.clear()
impaired.clear()
for w in range(len(y_all_healthy)):
    healthy.plot(y_all_healthy[w],color='blue')
for w in range(len(y_all_impaired)):
    impaired.plot(y_all_impaired[w],color='blue')
healthy.set_xlabel('n',fontsize=20)
healthy.set_ylabel('y',fontsize=20)
healthy.set_ylim(bottom=-20,top=140)
impaired.set_xlabel('n',fontsize=20)
impaired.set_ylabel('y',fontsize=20)
impaired.set_ylim(bottom=-20,top=140)

fig_healthy.savefig(path+'compare_healthy_impaired/all_clocks_overlaid/'+'/healthy_'+clock_type+'_clocks_y.png')
fig_impaired.savefig(path+'compare_healthy_impaired/all_clocks_overlaid/'+'/impaired_'+clock_type+'_clocks_y.png')

# just r
healthy.clear()
impaired.clear()
for w in range(len(r_all_healthy)):
    healthy.plot(r_all_healthy[w],color='blue')
for w in range(len(r_all_impaired)):
    impaired.plot(r_all_impaired[w],color='blue')
healthy.set_xlabel('n',fontsize=20)
healthy.set_ylabel('r',fontsize=20)
healthy.set_ylim(bottom=35,top=85)
impaired.set_xlabel('n',fontsize=20)
impaired.set_ylabel('r',fontsize=20)
impaired.set_ylim(bottom=35,top=85)

fig_healthy.savefig(path+'compare_healthy_impaired/all_clocks_overlaid/'+'/healthy_'+clock_type+'_clocks_r.png')
fig_impaired.savefig(path+'compare_healthy_impaired/all_clocks_overlaid/'+'/impaired_'+clock_type+'_clocks_r.png')

# just theta
healthy.clear()
impaired.clear()
for w in range(len(theta_all_healthy)):
    healthy.plot(theta_all_healthy[w],color='blue')
for w in range(len(theta_all_impaired)):
    impaired.plot(theta_all_impaired[w],color='blue')
healthy.set_xlabel('n',fontsize=20)
healthy.set_ylabel(r'$\theta$',fontsize=20)
healthy.set_ylim(bottom=0,top=6.5)
impaired.set_xlabel('n',fontsize=20)
impaired.set_ylabel(r'$\theta$',fontsize=20)
impaired.set_ylim(bottom=0,top=6.5)

fig_healthy.savefig(path+'compare_healthy_impaired/all_clocks_overlaid/'+'/healthy_'+clock_type+'_clocks_theta.png')
fig_impaired.savefig(path+'compare_healthy_impaired/all_clocks_overlaid/'+'/impaired_'+clock_type+'_clocks_theta.png')

# r vs. theta
healthy.clear()
impaired.clear()
for w in range(len(r_all_healthy)):
    healthy.plot(theta_all_healthy[w],r_all_healthy[w],color='blue')
for w in range(len(r_all_impaired)):
    impaired.plot(theta_all_impaired[w],r_all_impaired[w],color='blue')
healthy.set_xlabel(r'$\theta$',fontsize=20)
healthy.set_ylabel(r'$r$',fontsize=20)
healthy.set_xlim(left=0,right=6.5)
healthy.set_ylim(bottom=35,top=85)
impaired.set_xlabel(r'$\theta$',fontsize=20)
impaired.set_ylabel(r'$r$',fontsize=20)
impaired.set_xlim(left=0,right=6.5)
impaired.set_ylim(bottom=35,top=85)

fig_healthy.savefig(path+'compare_healthy_impaired/all_clocks_overlaid/'+'/healthy_'+clock_type+'_clocks_rtheta.png')
fig_impaired.savefig(path+'compare_healthy_impaired/all_clocks_overlaid/'+'/impaired_'+clock_type+'_clocks_rtheta.png')
