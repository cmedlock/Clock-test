# -*- coding: utf-8 -*-
# zoom in on the DFS coefficients at low frequencies

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

if not os.path.exists(path+'figs_raw'):
    os.makedirs(path+'figs_raw')

# save interesting quantities
energy_in_peak,energy_in_peak_command = [],[]

for fname in dirs[:4]:
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
    k = range(dft_size)
    dftx,dfty = np.fft.fft(x_zm,n=dft_size),np.fft.fft(y_zm,n=dft_size)
    
    # k_near_pi is the smallest k value for which w_k = 2*pi*k/N is
    # greater than pi
    k_near_pi,k_near_pi_command = 0,0
    if dft_size%2==0:
        k_near_pi = dft_size/2+1
    else:
        k_near_pi = math.ceil(dft_size/2)

    # only use the positive frequencies for plotting
    pos_k = k[:int(k_near_pi)]
    dftx_pos_k,dfty_pos_k = dftx[:k_near_pi],dfty[:k_near_pi]

    # zoom in on the DFS coefficients between ±(pi/8)

    # method 1: convert X[k] back to x[n] (can just use the original signal), 
    # extract the desired frequency band, then recalculate the DFS coefficients
        
    # convolve (truncated) sinc with x[n]
    omega_c = math.pi/8
    length_of_sinc = 601
    x1,y1 = [],[]
    for n in range(len(x_zm)+length_of_sinc-1):
        x1_of_n,y1_of_n = 0,0
        for m in range(len(x_zm)):
            x1_of_n += x_zm[m]*sinc(omega_c,n-(length_of_sinc-1)/2-m,length_of_sinc)
            y1_of_n += y_zm[m]*sinc(omega_c,n-(length_of_sinc-1)/2-m,length_of_sinc)
        x1.append(x1_of_n)
        y1.append(y1_of_n)
    # check dft of x1
    dftx1 = np.fft.fft(x1,n=len(x1))
    print 'check for x1: ',sum(x1),' = ',dftx1[0],' = ',8*dftx[0],'?'
    
    # decimate by 8 (d = decimated)
    x1_d,y1_d = [],[]
    for w in range(len(x1)):
        if w%8==0:
            x1_d.append(x1[w])
            y1_d.append(y1[w])
    
    x1_d,y1_d = np.array(x1_d),np.array(y1_d)
    
    # redo the N-point DFT
    dft_zoom_size = len(x1_d)
    dftx_zoom,dfty_zoom = np.fft.fft(x1_d,n=dft_zoom_size),np.fft.fft(y1_d,n=dft_zoom_size)
    #print 'check 1: ',sum(x1_d),' = ',dftx_zoom[0],' = ',dftx[0],'?'
    #print 'check 2: ',sum(y1_d),' = ',dfty_zoom[0],' = ',dfty[0],'?'
    #print 'check 3: ',sum(x1_d_command),' = ',dftx_zoom_command[0],' = ',dftx_command[0],'?'
    #print 'check 4: ',sum(y1_d_command),' = ',dfty_zoom_command[0],' = ',dfty_command[0],'?'

    # k_zoom_near_pi is the smallest k value for which w_k = 2*pi*k/N is
    # greater than pi
    k_zoom_near_pi = 0
    if dft_zoom_size%2==0:
        k_zoom_near_pi = dft_zoom_size/2+1
    else:
        k_zoom_near_pi = math.ceil(dft_zoom_size/2)

    # only use the positive frequencies for plotting
    pos_k_zoom = k_zoom[:int(k_zoom_near_pi)]
    dftx_zoom_pos_k,dfty_zoom_pos_k = dftx_zoom[:k_zoom_near_pi],dfty_zoom[:k_zoom_near_pi]

    # plot
    fig_xt,fig_yt = plt.figure(),plt.figure()
    fig_xt.subplots_adjust(hspace=0.5)
    fig_yt.subplots_adjust(hspace=0.5)
    xt,yt = fig_xt.add_subplot(311),fig_yt.add_subplot(311)
    xt.plot(t-t[0],x)
    yt.plot(t-t[0],y)
    
    dftxk,dftyk = fig_xt.add_subplot(312),fig_yt.add_subplot(312)    
    linecolor_dft = 'red'
    dftxk.plot(pos_k,np.abs(dftx_pos_k),linecolor_dft)
    dftyk.plot(pos_k,np.abs(dfty_pos_k),linecolor_dft)
    dftxk.set_ylim(bottom=min(np.abs(dftx_pos_k[1:]))/10,top=max(np.abs(dftx))*10)
    dftyk.set_ylim(bottom=min(np.abs(dfty_pos_k[1:]))/10,top=max(np.abs(dfty))*10)
    dftxk.set_yscale('log')
    dftyk.set_yscale('log')

    dftxk_zoom,dftyk_zoom = fig_xt.add_subplot(313),fig_yt.add_subplot(313)
    dftxk_zoom.plot(pos_k_zoom,np.abs(dftx_zoom_pos_k),linecolor_dft)
    dftyk_zoom.plot(pos_k_zoom,np.abs(dfty_zoom_pos_k),linecolor_dft)
    dftxk_zoom.set_ylim(bottom=min(np.abs(dftx_zoom_pos_k[1:]))/10,top=max(np.abs(dftx_zoom))*10)
    dftyk_zoom.set_ylim(bottom=min(np.abs(dfty_zoom_pos_k[1:]))/10,top=max(np.abs(dfty_zoom))*10)
    dftxk_zoom.set_yscale('log')
    dftyk_zoom.set_yscale('log')

    # set axis limits
    xt.set_xlim(right=max(t-t[0]))
    yt.set_xlim(right=max(t-t[0]))
    dftxk.set_xlim(right=pos_k[-1])
    dftyk.set_xlim(right=pos_k[-1])
    dftxk_zoom.set_xlim(right=pos_k_zoom[-1])
    dftyk_zoom.set_xlim(right=pos_k_zoom[-1])
    
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
    x_axis_fontsize = 15
    xt.set_xlabel('time [ms]',fontsize=x_axis_fontsize)
    yt.set_xlabel('time [ms]',fontsize=x_axis_fontsize)
    dftxk.set_xlabel(r'$\omega/(\frac{2\pi}{N})$',fontsize=x_axis_fontsize)
    dftyk.set_xlabel(r'$\omega/(\frac{2\pi}{N})$',fontsize=x_axis_fontsize)
    dftxk_zoom.set_xlabel(r'$8 \omega/(\frac{2\pi}{N})$',fontsize=x_axis_fontsize)
    dftyk_zoom.set_xlabel(r'$8 \omega/(\frac{2\pi}{N})$',fontsize=x_axis_fontsize)

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
