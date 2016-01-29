# 1.normalize for non-constant velocity
# 2.estimate ideal underlying sinusoids
# 3.compute correlation and mean squared difference 
# between ideal sinusoids and actual drawing
import math
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import numpy.fft
import os
from pylab import *

import clock_test as ct
ct = reload(ct)

path = '/Users/cmedlock/Documents/DSP_UROP/all_data/'
dirs = os.listdir(path)

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
    clock_type = 'COMMAND'
    x,y,t = [],[],[]

    # read in data
    record,found_clock = False,False
    for w in range(len(data)):
        line = data[w]
        # found copy clock?
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
                found_clock = False
                record = False
                continue
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
        for r in range(1,200):
            x_new = x[w]+r*dx/200
            y_new = y[w]+r*dy/200
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
            if abs(dist_total-dist_avg)<0.01:
                idx = j+1
	        break
        x_eqdist.append(x_interp[idx])
        y_eqdist.append(y_interp[idx])
    x_eqdist,y_eqdist = np.array(x_eqdist)-x_eqdist[0],np.array(y_eqdist)-y_eqdist[0]

    # now want to estimate the frequency of the underlying sinusoid
    # subtract mean values (zm = zero mean)
    x_eqdist_zm,y_eqdist_zm = x_eqdist-mean(x_eqdist),y_eqdist-mean(y_eqdist)
    
    # DFS coefficients
    dft_size = N_orig
    k = np.arange(dft_size)
    omega = 2.*k/float(dft_size)
    dftx,dfty = np.fft.fft(x_eqdist_zm,n=dft_size),np.fft.fft(y_eqdist,n=dft_size)

    # k_near_pi is the k value for which w_k = 2*pi*k/N is closest to,
    # but not larger than, pi
    k_near_pi = 0
    if dft_size%2==0:
        k_near_pi = dft_size/2+1
    else:
        k_near_pi = math.ceil(dft_size/2)

    # only use the positive frequencies
    posfreq = k[:k_near_pi]
    omega_posfreq = omega[:k_near_pi]
    dftx_posfreq,dfty_posfreq = dftx[:k_near_pi],dfty[:k_near_pi]
    #freq = 2*np.pi/dft_size*k
    #freq = np.fft.fftfreq(n=1024,d=1/(2*np.pi)) # NB: centered around 0
    
    # take the frequency of the largest DFS coefficient to be
    # the approximate frequency of the underlying sinusoid
    abs_dftx_posfreq,abs_dfty_posfreq = list(np.abs(dftx_posfreq)),list(np.abs(dfty_posfreq))
    k_true_x = abs_dftx_posfreq.index(max(abs_dftx_posfreq))
    k_true_y = abs_dfty_posfreq.index(max(abs_dfty_posfreq))
    w_true_x = 2*pi*k_true_x/dft_size
    w_true_y = 2*pi*k_true_y/dft_size
    #print '\nw_true_x = ',w_true_x,' and w_true_y = ',w_true_y
    #print '\n',np.angle(dftx[1]),np.angle(dftx[-1]),' ---> ',np.angle(dftx[1])-np.angle(dftx[-1])
    #print np.angle(dfty[1]),np.angle(dfty[-1]),' ---> ',np.angle(dfty[1])-np.angle(dfty[-1])

    # calculate ideal underlying sinusoid
    n = np.arange(len(x_eqdist))
    x_true,y_true = np.cos(w_true_x*n),np.sin(w_true_y*n)

    # use maximum correlation to determine phase and range of x or y to determine amplitude
    phase_x,max_corr_x = 0,0
    phase_y,max_corr_y = 0,0
    for w in range(len(x_true)):
        x_shifted = np.concatenate((x_true[w:],x_true[:w]))
        y_shifted = np.concatenate((y_true[w:],y_true[:w]))
        rho_x = np.dot(x_eqdist,x_shifted)
        rho_y = np.dot(y_eqdist,y_shifted)
        if rho_x>max_corr_x:
            max_corr_x = rho_x
            phase_x = w
        if rho_y>max_corr_y:
            max_corr_y = rho_y
            phase_y = w
    #print 'phase for x_true[n] is ',phase_x,' and max correlation is ',max_corr_x
    #print 'phase for y_true[n] is ',phase_y,' and max correlation is ',max_corr_y
    x_true = np.concatenate((x_true[phase_x:],x_true[:phase_x]))
    y_true = np.concatenate((y_true[phase_y:],y_true[:phase_y]))
    # use range of x or y to estimate amplitude
    amp_x = (max(x_eqdist)-min(x_eqdist))/2
    amp_y = (max(y_eqdist)-min(y_eqdist))/2
    x_true = amp_x*x_true
    y_true = amp_y*y_true
    x_true,y_true = x_true-x_true[0],y_true-y_true[0]
    
    # compute and save correlation between the two
    corr_x.append(np.dot(x_eqdist,x_true))
    corr_y.append(np.dot(y_eqdist,y_true))
    
    # compute and save mean squared difference between the two
    Etrue_x,Etrue_y = sum(x_true**2),sum(y_true**2)
    Ediff_x.append(sum((x_eqdist-x_true)**2))
    Ediff_y.append(sum((y_eqdist-y_true)**2))

    # figure settings
    linecolor_xy = 'blue'
    linecolor_dft = 'red'
    x_axis_fontsize = 15
    y_axis_fontsize = 20

    # initialize figures
    fig_xt,fig_yt,fig_xy = plt.figure(),plt.figure(),plt.figure()
    fig_xt.subplots_adjust(hspace=0.4,left=0.15)
    fig_yt.subplots_adjust(hspace=0.4,left=0.15)
    # declare subplots
    xt,yt,xy = fig_xt.add_subplot(311),fig_yt.add_subplot(311),fig_xy.add_subplot(111)
    xt_eqdist,yt_eqdist = fig_xt.add_subplot(312),fig_yt.add_subplot(312)
    dftxk,dftyk = fig_xt.add_subplot(313),fig_yt.add_subplot(313)
    # draw circles (actual and ideal)
    xy.plot(x,y,label='$x[n],y[n]$')
    xy.plot(x_true,y_true,'k-.',label='$x_{true}[n],y_{true}[n]$')
    xy.legend(loc='best',frameon=False)
    # draw x[n]
    xt.plot(x)
    xt_eqdist.plot(x_eqdist,label=r'$x_{eqdist}[n]$')
    xt_eqdist.plot(x_true,'k-.',lw=3,label=r'$x_{true}[n]$')
    xt_eqdist.legend(loc='best',frameon=False)
    dftxk.plot(omega_posfreq,np.abs(dftx_posfreq),linecolor_dft)
    # draw y[n]
    yt.plot(y)
    yt_eqdist.plot(y_eqdist,label=r'$y_{eqdist}[n]$')
    yt_eqdist.plot(y_true,'k-.',lw=3,label=r'$y_{true}[n]$')
    yt_eqdist.legend(loc='best',frameon=False)
    dftyk.plot(omega_posfreq,np.abs(dfty_posfreq),linecolor_dft)

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
    xt.set_xlim(right=len(x))
    xt_eqdist.set_xlim(right=len(x_eqdist))
    #dftxk.set_xlim(right=posfreq[-1])
    dftxk.set_ylim(bottom=min(np.abs(dftx_posfreq[1:]))/10,top=max(np.abs(dftx))*10)
    dftxk.set_yscale('log')

    yt.set_xlim(right=len(y))
    yt_eqdist.set_xlim(right=len(y_eqdist))
    #dftyk.set_xlim(right=posfreq[-1])
    dftyk.set_ylim(bottom=min(np.abs(dfty_posfreq[1:]))/10,top=max(np.abs(dfty))*10)
    dftyk.set_yscale('log')

    # set axis labels
    xy.set_xlabel(r'$x$',fontsize=x_axis_fontsize)
    xy.set_ylabel(r'$y$',fontsize=y_axis_fontsize)

    xt.set_xlabel(r'$n$',fontsize=x_axis_fontsize)
    xt.set_ylabel(r'$x[n]$',color=linecolor_xy,fontsize=y_axis_fontsize)
    xt_eqdist.set_xlabel(r'$n$',fontsize=x_axis_fontsize)
    xt_eqdist.set_ylabel('$x_{eqdist}[n]$,\n$x_{true}[n]$',fontsize=y_axis_fontsize-3)
    dftxk.set_xlabel(r'$\omega/\pi$',fontsize=x_axis_fontsize)
    for x1 in xt.get_yticklabels():
        x1.set_color(linecolor_xy)
    for x2 in xt_eqdist.get_yticklabels():
        x2.set_color(linecolor_xy)

    yt.set_xlabel(r'$n$',fontsize=x_axis_fontsize)
    yt.set_ylabel(r'$y[n]$',color=linecolor_xy,fontsize=y_axis_fontsize)
    yt_eqdist.set_xlabel(r'$n$',fontsize=x_axis_fontsize)
    yt_eqdist.set_ylabel('$y_{eqdist}[n]$,\n$y_{true}[n]$',fontsize=y_axis_fontsize-3)
    dftyk.set_xlabel(r'$\omega/\pi$',fontsize=x_axis_fontsize)
    for y1 in yt.get_yticklabels():
        y1.set_color(linecolor_xy)
    for y2 in yt_eqdist.get_yticklabels():
        y2.set_color(linecolor_xy)
        
    dftxk.set_ylabel(r'$|X_{eqdist}[k]|$',color=linecolor_dft,fontsize=y_axis_fontsize)
    dftyk.set_ylabel(r'$|Y_{eqdist}[k]|$',color=linecolor_dft,fontsize=y_axis_fontsize)
    for v1 in dftxk.get_yticklabels():
        v1.set_color(linecolor_dft)
    for v2 in dftyk.get_yticklabels():
        v2.set_color(linecolor_dft)
    
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
    fig_xt.savefig(path+'figs_raw/'+fname[:len(fname)-4]+'/x_true_'+clock_type+'_'+fname[:len(fname)-4]+'.png')
    fig_yt.savefig(path+'figs_raw/'+fname[:len(fname)-4]+'/y_true_'+clock_type+'_'+fname[:len(fname)-4]+'.png')

    plt.close('all')
