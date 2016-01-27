# -*- coding: utf-8 -*-
# for sanity checks on zooming in on the low frequency DFS coefficients

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

for fname in dirs[:3]:
    
    if 'Scored' not in fname:
        continue
    
    if not os.path.exists(path+'figs_raw/'+fname[:len(fname)-4]):
        os.makedirs(path+'figs_raw/'+fname[:len(fname)-4])

    # get coordinates and timestamps
    n = np.arange(250)
    x = np.cos(2*np.pi/250*n)
    y = np.sin(2*np.pi/250*n)
        
    # plot spectrum of x[t] and y[t]
    # subtract mean values (zm = zero mean)
    x_zm,y_zm = x-mean(x),y-mean(y)
    # DFS coefficients
    dft_size = len(x_zm)
    dftx,dfty = np.fft.fft(x_zm,n=dft_size),np.fft.fft(y_zm,n=dft_size)
    k = np.arange(dft_size)
    freq = 2*np.pi/dft_size*k
    
    # k_near_pi is the smallest k value for which w_k = 2*pi*k/N is
    # greater than pi
    k_near_pi = 0
    if dft_size%2==0:
        k_near_pi = dft_size/2+1
    else:
        k_near_pi = math.ceil(dft_size/2)

    # only use the positive frequencies for plotting
    pos_k = k[:int(k_near_pi)]
    dftx_pos_k,dfty_pos_k = dftx[:k_near_pi],dfty[:k_near_pi]

    # zoom in on the DFS coefficients between ±(pi/8)

    # convolve (truncated) sinc with x[n]
    omega_c = math.pi/8
    length_of_sinc = 601
    x1,y1 = [],[]
    for w in range(len(x_zm)+length_of_sinc-1):
        x1_of_w,y1_of_w = 0,0
        for m in range(len(x_zm)):
            x1_of_w += x_zm[m]*sinc(omega_c,w-(length_of_sinc-1)/2-m,length_of_sinc)
            y1_of_w += y_zm[m]*sinc(omega_c,w-(length_of_sinc-1)/2-m,length_of_sinc)
        x1.append(x1_of_w)
        y1.append(y1_of_w)
    # check filter
    h1 = []
    for w in range(-(length_of_sinc-1)/2,(length_of_sinc-1)/2+1):
        h1.append(sinc(omega_c,w,length_of_sinc))
    dfth1 = np.fft.fft(h1,n=dft_size)
    kh1 = np.arange(dft_size)
    # check dft of x1
    dftx1 = np.fft.fft(x1,n=dft_size)
    k1 = np.arange(dft_size)
    print 'check for x1: ',sum(x1),' = ',np.abs(dftx1[0]),' = ',np.abs(8*dftx[0]),'?'
    print np.abs(dftx)[:3]
    print np.abs(dftx1)[:3]
    print np.abs(dftx[35]*dfth1[35]),np.abs(dftx1[35])
    
    # decimate by 8 (d = decimated)
    x1_d,y1_d = [],[]
    for w in range(len(x1)):
        if w%8==0:
            x1_d.append(x1[w])
            y1_d.append(y1[w])

    x1_d,y1_d = np.array(x1_d),np.array(y1_d)

    # redo the N-point DFT
    dft_zoom_size = len(x1_d)
    k_zoom = range(dft_zoom_size)
    dftx_zoom,dfty_zoom = np.fft.fft(x1_d,n=dft_zoom_size),np.fft.fft(y1_d,n=dft_zoom_size)
    print 'check 1: ',sum(x1_d),' = ',np.abs(dftx_zoom[0]),' = ',np.abs(dftx[0]),'?'
    #print np.abs(dftx_zoom)[:10]
    #print 'check 2: ',np.abs(max(dftx_zoom)),' = ',np.abs(max(dftx)),'?'

    # test plot for zoom (don't save)
    plt.close('all')
    fig = plt.figure()
    ax1 = fig.add_subplot(411)
    ax1.stem(k,np.abs(dftx))
    ax1.set_ylabel('|X[k]|')
    ax1.set_xlim(left=-1,right=max(k)+1)
    ax1.set_ylim(bottom=-1,top=max(abs(dftx))+1)
    ax2 = fig.add_subplot(412)
    ax2.stem(kh1,np.abs(dfth1))
    ax2.set_ylabel('|H_1[k]|')
    ax2.set_xlim(left=-1,right=max(kh1)+1)
    ax2.set_ylim(bottom=-1,top=max(np.abs(dfth1))+1)
    ax3 = fig.add_subplot(413)
    ax3.stem(k1,np.abs(dftx1))
    ax3.set_ylabel('|X_1[k]|')
    ax3.set_xlim(left=-1,right=max(k1)+1)
    ax3.set_ylim(bottom=-1,top=max(abs(dftx1))+1)
    ax4 = fig.add_subplot(414)
    ax4.stem(k_zoom,np.abs(dftx_zoom))
    ax4.set_ylabel('|X_zoom[k]|')
    ax4.set_xlim(left=-1,right=max(k_zoom)+1)
    ax4.set_ylim(bottom=-1,top=max(np.abs(dftx_zoom))+1)
    plt.show()
"""
    # k_zoom_near_pi is the smallest k value for which w_k = 2*pi*k/N is
    # greater than or equal to pi
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
    xt.plot(n,x)
    yt.plot(n,y)
    
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
    xt.set_xlim(right=max(n))
    yt.set_xlim(right=max(n))
    dftxk.set_xlim(right=pos_k[-1])
    dftyk.set_xlim(right=pos_k[-1])
    dftxk_zoom.set_xlim(right=pos_k_zoom[-1])
    dftyk_zoom.set_xlim(right=pos_k_zoom[-1])
    
    # set titles and file labels
    title_fontsize = 20
    #xt.set_title('Sample Copy Clock',fontsize=title_fontsize)
    #yt.set_title('Sample Copy Clock',fontsize=title_fontsize)
    #xy.set_title('Sample Copy Clock',fontsize=title_fontsize)

    # set axis labels
    x_axis_fontsize = 15
    xt.set_xlabel('n',fontsize=x_axis_fontsize)
    yt.set_xlabel('n',fontsize=x_axis_fontsize)
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
    fig_xt.savefig(path+'figs_raw/healthy_robot/dftx_healthy_robot.png')
    fig_yt.savefig(path+'figs_raw/healthy_robot/dfty_healthy_robot.png')    
"""
print 'Done'
