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

path = '/Users/cmedlock/Documents/DSP_UROP/all_data/'
dirs = os.listdir(path)

if not os.path.exists(path+'figs_raw'):
    os.makedirs(path+'figs_raw')

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

    x,y,t = np.array(x),np.array(y),np.array(t)

    # subtract mean values (zm = zero mean)
    x_zm,y_zm = x-mean(x),y-mean(y)
    
    # DFS coefficients
    dft_size = len(x)
    dftx,dfty = np.fft.fft(x_zm,n=dft_size),np.fft.fft(y_zm,n=dft_size)
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
    
    # copy clocks
    fig_xt,fig_yt = plt.figure(),plt.figure()
    fig_xt.subplots_adjust(hspace=0.3)
    fig_yt.subplots_adjust(hspace=0.3)
    xt,yt = fig_xt.add_subplot(211),fig_yt.add_subplot(211)
    xt.plot(t-t[0],x)
    yt.plot(t-t[0],y)
    
    dftxk,dftyk = fig_xt.add_subplot(212),fig_yt.add_subplot(212)    
    linecolor_dft = 'red'
    dftxk.plot(posfreq,np.abs(dftx_posfreq),linecolor_dft)
    dftyk.plot(posfreq,np.abs(dfty_posfreq),linecolor_dft)
    dftxk.set_ylim(bottom=min(np.abs(dftx_posfreq[1:]))/10,top=max(np.abs(dftx))*10)
    dftyk.set_ylim(bottom=min(np.abs(dfty_posfreq[1:]))/10,top=max(np.abs(dfty))*10)
    dftxk.set_yscale('log')
    dftyk.set_yscale('log')

    # set axis limits
    xt.set_xlim(right=max(t-t[0]))
    yt.set_xlim(right=max(t-t[0]))
    dftxk.set_xlim(right=posfreq[-1])
    dftyk.set_xlim(right=posfreq[-1])
    
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
    xt.set_xlabel('time [ms]',fontsize=x_axis_fontsize)
    yt.set_xlabel('time [ms]',fontsize=x_axis_fontsize)
    dftxk.set_xlabel(r'$k$',fontsize=x_axis_fontsize)
    dftyk.set_xlabel(r'$k$',fontsize=x_axis_fontsize)

    y_axis_fontsize = 25
    linecolor_xy = 'blue'
    xt.set_ylabel(r'$x$',color=linecolor_xy,fontsize=y_axis_fontsize)
    yt.set_ylabel(r'$y$',color=linecolor_xy,fontsize=y_axis_fontsize)
    for x1 in xt.get_yticklabels():
        x1.set_color(linecolor_xy)
    for y1 in yt.get_yticklabels():
        y1.set_color(linecolor_xy)
        
    dftxk.set_ylabel(r'$|X[k]|$',color=linecolor_dft,fontsize=y_axis_fontsize)
    dftyk.set_ylabel(r'$|Y[k]|$',color=linecolor_dft,fontsize=y_axis_fontsize)
    for v1 in dftxk.get_yticklabels():
        v1.set_color(linecolor_dft)
    for v2 in dftyk.get_yticklabels():
        v2.set_color(linecolor_dft)
    
    # save figures
    fig_xt.savefig(path+'figs_raw/'+fname[:len(fname)-4]+'/dftx_'+clock_type+'_'+fname[:len(fname)-4]+'.png')
    fig_yt.savefig(path+'figs_raw/'+fname[:len(fname)-4]+'/dfty_'+clock_type+'_'+fname[:len(fname)-4]+'.png')

    plt.close('all')
