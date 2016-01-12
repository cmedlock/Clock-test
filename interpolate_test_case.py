# -*- coding: utf-8 -*-
import math
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import numpy.fft
import os
from pylab import *

# length of truncated sinc (longer = more accurate interpolation)
length = 100
missing_points = [5,15,27]

# original signal, appropriately bandlimited so that
# it can be downsampled by 2 with no aliasing,
# i.e. X(e^jw)=0 for |w|>= pi/2
n = np.arange(length)
x = np.sin(np.pi/3*n)
# remove some points
for idx in missing_points:
    x[idx] = -5
# interpolate to find the missing points
for point in missing_points:
    # downsample then expand by 2
    # if the missing point had an even index, want to save all of the odd samples
    if point%2==0:
        x[::2] = 0
        expanded = x
    # if the missing point had an odd index, want to save all of the even samples
    elif point%2==1:
        pass
"""
    # interpolate to find missing points if necessary
    # assume for right now that we are never missing multiple consecutive points
    sinc_copy = None
    xexpand_copy,yexpand_copy = [],[]
    xexpand_copy_forplot,yexpand_copy_forplot = [],[] # for the expanded signal with no interpolation
    if -5 in x_copy:
        # expand by 2 (insert zeros)
        print '   COPY clock has missing points, expanding x[t] and y[t]...'
        for v in range(len(x_copy)-1):
            xexpand_copy.append(x_copy[v])
            yexpand_copy.append(y_copy[v])
            if x_copy[v]!=-5 and x_copy[v+1]!=-5:
                xexpand_copy.append(0)
                yexpand_copy.append(0)
        xexpand_copy.append(x_copy[-1])
        yexpand_copy.append(y_copy[-1])
        xexpand_copy_forplot = [elt for elt in xexpand_copy]
        yexpand_copy_forplot = [elt for elt in yexpand_copy]
        for d in range(len(xexpand_copy)):
            if xexpand_copy[d]!=-5:
                continue
            print '      found missing point at index ',(d-1)/2
            # windowed sinc filter
            print '         forming windowed sinc filter...'
            L = 2
            n_neg = np.arange(-d,0)
            n_pos = np.arange(1,len(xexpand_copy)-d)
            sinc_neg = np.sin(np.pi*n_neg/L)/(np.pi*n_neg/L)
            sinc_pos = np.sin(np.pi*n_pos/L)/(np.pi*n_pos/L)
            # NB: normally would set sinc[0] equal to 1, but in this
            # case we don't want the -5 at that point to contribute to
            # the sum, so just set sinc[0] equal to 0
            sinc_copy = np.concatenate((sinc_neg,np.array([0]),sinc_pos))
            # evaluate convolution at missing point
            missing_x = np.dot(xexpand_copy,sinc_copy)
            missing_y = np.dot(yexpand_copy,sinc_copy)
            print '         missing x is ',missing_x,' and missing y is ',missing_y
            xexpand_copy[d] = missing_x
            yexpand_copy[d] = missing_y
        # downsample by 2 to remove the zeros
        xintp_copy,yintp_copy = [],[]
        for v in range(len(xexpand_copy)):
            if v%2==0:
                xintp_copy.append(xexpand_copy[v])
                yintp_copy.append(yexpand_copy[v])
            elif v%2==1 and xexpand_copy[v]!=0:
                print '      inserting missing point (',xexpand_copy[v],',',yexpand_copy[v],') at index ',v/2
                xintp_copy.append(xexpand_copy[v])
                yintp_copy.append(yexpand_copy[v])
        print len(x_copy),len(y_copy)
        print len(xintp_copy),len(yintp_copy)

    xexpand_command,yexpand_command = [],[]
    xexpand_command_forplot,yexpand_command_forplot = [],[] # for the expanded signal with no interpolation
    sinc_command = 0
    if -5 in x_command:
        # expand by 2 (insert zeros)
        print '   COMMAND clock has missing points, expanding x[t] and y[t]...'
        for v in range(len(x_command)-1):
            xexpand_command.append(x_command[v])
            yexpand_command.append(y_command[v])
            if x_command[v]!=-5 and x_command[v+1]!=-5:
                xexpand_command.append(0)
                yexpand_command.append(0)
        xexpand_command.append(x_command[-1])
        yexpand_command.append(y_command[-1])
        xexpand_command_forplot = [elt for elt in xexpand_command]
        yexpand_command_forplot = [elt for elt in yexpand_command]
        for d in range(len(xexpand_command)):
            if xexpand_command[d]!=-5:
                continue
            print '      found missing point at index ',(d-1)/2
            # windowed sinc filter
            print '         forming windowed sinc filter...'
            L = 2
            n_neg = np.arange(-d,0)
            n_pos = np.arange(1,len(xexpand_command)-d)
            sinc_neg = np.sin(np.pi*n_neg/L)/(np.pi*n_neg/L)
            sinc_pos = np.sin(np.pi*n_pos/L)/(np.pi*n_pos/L)
            # NB: normally would set sinc[0] equal to 1, but in this
            # case we don't want the -5 at that point to contribute to
            # the sum, so just set sinc[0] equal to 0
            sinc_command = np.concatenate((sinc_neg,np.array([0]),sinc_pos))
            # evaluate convolution at missing point
            missing_x = np.dot(xexpand_command,sinc_command)
            missing_y = np.dot(yexpand_command,sinc_command)
            print '         missing x is ',missing_x,' and missing y is ',missing_y
            xexpand_command[d] = missing_x
            yexpand_command[d] = missing_y
        # downsample by 2 to remove the zeros
        xintp_command,yintp_command = [],[]
        for v in range(len(xexpand_command)):
            if v%2==0:
                xintp_command.append(xexpand_command[v])
                yintp_command.append(yexpand_command[v])
            elif v%2==1 and xexpand_command[v]!=0:
                print '      inserting missing point (',xexpand_command[v],',',yexpand_command[v],') at index ',v/2
                xintp_command.append(xexpand_command[v])
                yintp_command.append(yexpand_command[v])
        print len(x_command),len(y_command)
        print len(xintp_command),len(yintp_command)

    # copy clocks
    fig_xt_copy,fig_yt_copy = plt.figure(),plt.figure()
    # no interpolation
    xt_copy,yt_copy = fig_xt_copy.add_subplot(311),fig_yt_copy.add_subplot(311)
    xt_copy.stem(np.arange(30),xexpand_copy_forplot[:30],linestyle='none',marker='o')
    yt_copy.stem(np.arange(30),yexpand_copy_forplot[:30],linestyle='none',marker='o')
    xt_copy.set_ylim(bottom=-10,top=max(xintp_copy[:30])+10)
    yt_copy.set_ylim(bottom=-10,top=max(yintp_copy[:30])+10)
    # windowed sinc filter for the last missing point
    sincxt_copy,sincyt_copy = fig_xt_copy.add_subplot(312),fig_yt_copy.add_subplot(312)
    sincxt_copy.stem(np.arange(30),sinc_copy[:30],linestyle='none',marker='o')
    sincyt_copy.stem(np.arange(30),sinc_copy[:30],linestyle='none',marker='o')
    # interpolated signal
    xintpt_copy,yintpt_copy = fig_xt_copy.add_subplot(313),fig_yt_copy.add_subplot(313)
    xintpt_copy.stem(np.arange(30),xexpand_copy[:30],linestyle='none',marker='o')
    yintpt_copy.stem(np.arange(30),yexpand_copy[:30],linestyle='none',marker='o')
    xintpt_copy.set_ylim(bottom=-10,top=max(xintp_copy[:30])+10)
    yintpt_copy.set_ylim(bottom=-10,top=max(yintp_copy[:30])+10)
    
    # set titles and file labels
    title_fontsize = 20
    #xt_copy.set_title('Sample Copy Clock',fontsize=title_fontsize)
    #yt_copy.set_title('Sample Copy Clock',fontsize=title_fontsize)
    #xy_copy.set_title('Sample Copy Clock',fontsize=title_fontsize)

    fig_xt_copy.text(0.99, 0.96, fname[:len(fname)-4],fontsize=10,color='red',va='baseline',ha='right',multialignment='left')
    fig_yt_copy.text(0.99, 0.96, fname[:len(fname)-4],fontsize=10,color='red',va='baseline',ha='right',multialignment='left')

    if 'YDU' in fname:
        fig_xt_copy.text(0.25, 0.955, 'HEALTHY',fontsize=15,color='black',va='baseline',ha='right',multialignment='left')
        fig_yt_copy.text(0.25, 0.955, 'HEALTHY',fontsize=15,color='black',va='baseline',ha='right',multialignment='left')
    elif 'CIN' in fname:
        fig_xt_copy.text(0.25, 0.955, 'IMPAIRED',fontsize=15,color='black',va='baseline',ha='right',multialignment='left')
        fig_yt_copy.text(0.25, 0.955, 'IMPAIRED',fontsize=15,color='black',va='baseline',ha='right',multialignment='left')
    else:
        print 'not a valid filename'

    # set axis labels
    x_axis_fontsize = 20
    xintpt_copy.set_xlabel('n',fontsize=x_axis_fontsize)
    yintpt_copy.set_xlabel('n',fontsize=x_axis_fontsize)

    y_axis_fontsize = 20
    xt_copy.set_ylabel('x[n] \n (w/out interpolation)',fontsize=y_axis_fontsize-7)
    yt_copy.set_ylabel('y[n] \n (w/out interpolation)',fontsize=y_axis_fontsize-7)
    sincxt_copy.set_ylabel('sinc[n]',fontsize=y_axis_fontsize)
    sincyt_copy.set_ylabel('sinc[n]',fontsize=y_axis_fontsize)
    xintpt_copy.set_ylabel('x[n] \n (w/ interpolation)',fontsize=y_axis_fontsize-7)
    yintpt_copy.set_ylabel('y[n] \n (w/ interpolation)',fontsize=y_axis_fontsize-7)

    # save figures
    fig_xt_copy.savefig(path+'figs_raw/'+fname[:len(fname)-4]+'/xintpt_copy_'+fname[:len(fname)-4]+'.png')
    fig_yt_copy.savefig(path+'figs_raw/'+fname[:len(fname)-4]+'/yintpt_copy_'+fname[:len(fname)-4]+'.png')

    # command clocks
    fig_xt_command,fig_yt_command = plt.figure(),plt.figure()
    # no interpolation
    xt_command,yt_command = fig_xt_command.add_subplot(311),fig_yt_command.add_subplot(311)
    xt_command.stem(np.arange(30),xexpand_command_forplot[:30],linestyle='none',marker='o')
    yt_command.stem(np.arange(30),yexpand_command_forplot[:30],linestyle='none',marker='o')
    xt_command.set_ylim(bottom=-10,top=max(xintp_command[:30])+10)
    yt_command.set_ylim(bottom=-10,top=max(yintp_command[:30])+10)
    # windowed sinc filter for the last missing point
    sincxt_command,sincyt_command = fig_xt_command.add_subplot(312),fig_yt_command.add_subplot(312)
    sincxt_command.stem(np.arange(30),sinc_command[:30],linestyle='none',marker='o')
    sincyt_command.stem(np.arange(30),sinc_command[:30],linestyle='none',marker='o')
    # interpolated signal
    xintpt_command,yintpt_command = fig_xt_command.add_subplot(313),fig_yt_command.add_subplot(313)
    xintpt_command.stem(np.arange(30),xexpand_command[:30],linestyle='none',marker='o')
    yintpt_command.stem(np.arange(30),yexpand_command[:30],linestyle='none',marker='o')
    xintpt_command.set_ylim(bottom=-10,top=max(xintp_command[:30])+10)
    yintpt_command.set_ylim(bottom=-10,top=max(yintp_command[:30])+10)
    
    # set titles and file labels
    title_fontsize = 20
    #xt_command.set_title('Sample command Clock',fontsize=title_fontsize)
    #yt_command.set_title('Sample command Clock',fontsize=title_fontsize)
    #xy_command.set_title('Sample command Clock',fontsize=title_fontsize)

    fig_xt_command.text(0.99, 0.96, fname[:len(fname)-4],fontsize=10,color='red',va='baseline',ha='right',multialignment='left')
    fig_yt_command.text(0.99, 0.96, fname[:len(fname)-4],fontsize=10,color='red',va='baseline',ha='right',multialignment='left')

    if 'YDU' in fname:
        fig_xt_command.text(0.25, 0.955, 'HEALTHY',fontsize=15,color='black',va='baseline',ha='right',multialignment='left')
        fig_yt_command.text(0.25, 0.955, 'HEALTHY',fontsize=15,color='black',va='baseline',ha='right',multialignment='left')
    elif 'CIN' in fname:
        fig_xt_command.text(0.25, 0.955, 'IMPAIRED',fontsize=15,color='black',va='baseline',ha='right',multialignment='left')
        fig_yt_command.text(0.25, 0.955, 'IMPAIRED',fontsize=15,color='black',va='baseline',ha='right',multialignment='left')
    else:
        print 'not a valid filename'

    # set axis labels
    x_axis_fontsize = 20
    xintpt_command.set_xlabel('n',fontsize=x_axis_fontsize)
    yintpt_command.set_xlabel('n',fontsize=x_axis_fontsize)

    y_axis_fontsize = 20
    xt_command.set_ylabel('x[n] \n (w/out interpolation)',fontsize=y_axis_fontsize-7)
    yt_command.set_ylabel('y[n] \n (w/out interpolation)',fontsize=y_axis_fontsize-7)
    sincxt_command.set_ylabel('sinc[n]',fontsize=y_axis_fontsize)
    sincyt_command.set_ylabel('sinc[n]',fontsize=y_axis_fontsize)
    xintpt_command.set_ylabel('x[n] \n (w/ interpolation)',fontsize=y_axis_fontsize-7)
    yintpt_command.set_ylabel('y[n] \n (w/ interpolation)',fontsize=y_axis_fontsize-7)

    # save figures
    fig_xt_command.savefig(path+'figs_raw/'+fname[:len(fname)-4]+'/xintpt_command_'+fname[:len(fname)-4]+'.png')
    fig_yt_command.savefig(path+'figs_raw/'+fname[:len(fname)-4]+'/yintpt_command_'+fname[:len(fname)-4]+'.png')

    plt.close('all')
"""