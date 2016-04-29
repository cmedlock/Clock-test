# plot DFS coefficients of each x vs. t and y vs. t signal

import math
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import numpy.fft
import os
from pylab import *

import clock_test as ct
ct = reload(ct)

path = '/Users/cmedlock/Documents/DSP_UROP/DataCatherine/'
dirs = os.listdir(path)

if not os.path.exists(path+'figs_raw'):
    os.makedirs(path+'figs_raw')

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
    if len(x)==0:
        continue
        
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
    
    # DFS coefficients
    dft_size = len(x_eqdist)
    dftx,dfty = np.fft.fft(x_eqdist,n=dft_size),np.fft.fft(y_eqdist,n=dft_size)
    k = np.arange(dft_size)

    # k_near_pi is the k value for which w_k = 2*pi*k/N is closest to,
    # but not larger than, pi
    k_near_pi = 0
    if dft_size%2==0:
        k_near_pi = dft_size/2+1
    else:
        k_near_pi = math.ceil(dft_size/2)

    # only use the positive frequencies
    posfreq = k[:k_near_pi]
    dftx_posfreq,dfty_posfreq = dftx[:k_near_pi],dfty[:k_near_pi]
    #freq = 2*np.pi/dft_size*k
    #freq = np.fft.fftfreq(n=1024,d=1/(2*np.pi)) # NB: centered around 0
    
    # plot
    fig_xt,fig_yt = plt.figure(),plt.figure()
    fig_xt.subplots_adjust(hspace=0.3)
    fig_yt.subplots_adjust(hspace=0.3)
    
    dftxk,dftyk = fig_xt.add_subplot(211),fig_yt.add_subplot(211)    
    linecolor_dft = 'red'
    dftxk.plot(posfreq,np.abs(dftx_posfreq),linecolor_dft)
    dftyk.plot(posfreq,np.abs(dfty_posfreq),linecolor_dft)
    dftxk.set_ylim(bottom=min(np.abs(dftx_posfreq[1:]))/10,top=max(np.abs(dftx))*10)
    dftyk.set_ylim(bottom=min(np.abs(dfty_posfreq[1:]))/10,top=max(np.abs(dfty))*10)
    dftxk.set_yscale('log')
    dftyk.set_yscale('log')

    dftxk_zoom,dftyk_zoom = fig_xt.add_subplot(212),fig_yt.add_subplot(212)
    dftxk_zoom.stem(posfreq[:10],np.abs(dftx_posfreq[:10]),linecolor_dft,markerfmt='ro')
    dftyk_zoom.stem(posfreq[:10],np.abs(dfty_posfreq[:10]),linecolor_dft,markerfmt='ro')
    dftxk_zoom.set_ylim(bottom=min(np.abs(dftx_posfreq[1:10]))/10,top=max(np.abs(dftx))*10)
    dftyk_zoom.set_ylim(bottom=min(np.abs(dfty_posfreq[1:10]))/10,top=max(np.abs(dfty))*10)
    dftxk_zoom.set_yscale('log')
    dftyk_zoom.set_yscale('log')
    
    # set axis limits
    dftxk.set_xlim(right=posfreq[-1])
    dftyk.set_xlim(right=posfreq[-1])
    dftxk_zoom.set_xlim(left=0,right=posfreq[10])
    dftyk_zoom.set_xlim(left=0,right=posfreq[10])
    
    # set titles and file labels
    title_fontsize = 20
    #xt.set_title('Sample Copy Clock',fontsize=title_fontsize)
    #yt.set_title('Sample Copy Clock',fontsize=title_fontsize)
    #xy.set_title('Sample Copy Clock',fontsize=title_fontsize)

    fig_xt.text(0.99, 0.96,fname[:len(fname)-4],fontsize=10,color='red',va='baseline',ha='right',multialignment='left')
    fig_yt.text(0.99, 0.96,fname[:len(fname)-4],fontsize=10,color='red',va='baseline',ha='right',multialignment='left')

    if 'YDU' in fname:
        fig_xt.text(0.25, 0.955, 'HEALTHY',fontsize=15,color='black',va='baseline',ha='right',multialignment='left')
        fig_yt.text(0.25, 0.955, 'HEALTHY',fontsize=15,color='black',va='baseline',ha='right',multialignment='left')
    elif 'CIN' in fname:
        fig_xt.text(0.25, 0.955, 'IMPAIRED',fontsize=15,color='black',va='baseline',ha='right',multialignment='left')
        fig_yt.text(0.25, 0.955, 'IMPAIRED',fontsize=15,color='black',va='baseline',ha='right',multialignment='left')
    else:
        print 'not a valid filename'
        
    # set axis labels
    x_axis_fontsize = 20
    #xt.set_xlabel('time [ms]',fontsize=x_axis_fontsize)
    #yt.set_xlabel('time [ms]',fontsize=x_axis_fontsize)
    dftxk.set_xlabel(r'$k$',fontsize=x_axis_fontsize)
    dftyk.set_xlabel(r'$k$',fontsize=x_axis_fontsize)
    dftxk_zoom.set_xlabel(r'$k$',fontsize=x_axis_fontsize)
    dftyk_zoom.set_xlabel(r'$k$',fontsize=x_axis_fontsize)

    y_axis_fontsize = 25
    linecolor_xy = 'blue'
    #xt.set_ylabel(r'$x$',color=linecolor_xy,fontsize=y_axis_fontsize)
    #yt.set_ylabel(r'$y$',color=linecolor_xy,fontsize=y_axis_fontsize)
    #for x1 in xt.get_yticklabels():
    #    x1.set_color(linecolor_xy)
    #for y1 in yt.get_yticklabels():
    #    y1.set_color(linecolor_xy)
        
    dftxk.set_ylabel(r'$|X[k]|$',color=linecolor_dft,fontsize=y_axis_fontsize)
    dftyk.set_ylabel(r'$|Y[k]|$',color=linecolor_dft,fontsize=y_axis_fontsize)
    for v1 in dftxk.get_yticklabels():
        v1.set_color(linecolor_dft)
    for v2 in dftyk.get_yticklabels():
        v2.set_color(linecolor_dft)
    dftxk_zoom.set_ylabel(r'$|X[k]|$',color=linecolor_dft,fontsize=y_axis_fontsize)
    dftyk_zoom.set_ylabel(r'$|Y[k]|$',color=linecolor_dft,fontsize=y_axis_fontsize)
    for v1 in dftxk_zoom.get_yticklabels():
        v1.set_color(linecolor_dft)
    for v2 in dftyk_zoom.get_yticklabels():
        v2.set_color(linecolor_dft)
    
    # save figures
    fig_xt.savefig(path+'figs_raw/'+fname[:len(fname)-4]+'/dftx_'+clock_type+'_'+fname[:len(fname)-4]+'.png')
    fig_yt.savefig(path+'figs_raw/'+fname[:len(fname)-4]+'/dfty_'+clock_type+'_'+fname[:len(fname)-4]+'.png')

    plt.close('all')
