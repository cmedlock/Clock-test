# -*- coding: utf-8 -*-
import math
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import numpy.fft
import os
from pylab import *

path = '/Users/cmedlock/Documents/DSP_UROP/all_data/'
dirs = os.listdir(path)

if not os.path.exists(path+'figs_raw'):
    os.makedirs(path+'figs_raw')

for fname in dirs[:3]:
    if 'Scored' not in fname:
        continue
    print 'reading file ',fname,'...'
    f = open(path+fname)
    data = f.readlines()

    if not os.path.exists(path+'figs_raw/'+fname[:len(fname)-4]):
        os.makedirs(path+'figs_raw/'+fname[:len(fname)-4])

    # copy or command clock?
    clock_type = ''
    x_copy,y_copy,p_copy,t_copy = [],[],[],[]
    x_command,y_command,p_command,t_command = [],[],[],[]

    # read in data
    record,found_clock = False,False
    for w in range(len(data)):
        line = data[w]
        # found copy clock?
        if found_clock==False:
            if 'COPY' in line:
                clock_type = 'COPY'
                found_clock = True
            elif 'COMMAND' in line:
                clock_type = 'COMMAND'
                found_clock = True
            continue
        # start recording?
        elif found_clock==True and 'CLOCKFACE' in line:
            record = True
            continue
        # stop recording?
        elif record==True:
            if 'symbol label' in line and clock_type=='COPY' and len(x_copy)>0:
                found_clock = False
                record = False
                continue
            elif 'symbol label' in line and clock_type=='COMMAND' and len(x_command)>0:
                found_clock = False
                record = False
                continue
            elif 'point' not in line:
                continue
        # other?
        else:
            continue
        line = line.split('"')
        
        # position and time
        xcoord = double(line[3])
        ycoord = double(line[1])
        pressure = double(line[5])
        timestamp = double(line[7])
        if clock_type=='COPY':
            # if a point has been skipped, leave a space for it with the appropriate timestamp
            if len(x_copy)>0 and timestamp-t_copy[-1]>20:
                x_copy.append(-5)
                y_copy.append(-5)
                p_copy.append(-15)
                t_copy.append(t_copy[-1]+(timestamp-t_copy[-1])/2)
            x_copy.append(xcoord)
            y_copy.append(ycoord)
            p_copy.append(pressure)
            t_copy.append(timestamp)
        elif clock_type=='COMMAND':
            if len(x_command)>0 and timestamp-t_command[-1]>20:
                x_command.append(-5)
                y_command.append(-5)
                p_command.append(-5)
                t_command.append(t_command[-1]+(timestamp-t_command[-1])/2)
            x_command.append(xcoord)
            y_command.append(ycoord)
            p_command.append(pressure)
            t_command.append(timestamp)
        else:
            print 'not a valid clock type'
        
        try:
            if t_copy[-1]-t_copy[-2]>19:
                print '   ',w,found_clock,clock_type,record,xcoord,ycoord,t_copy[-1]-t_copy[-2]
        except:
            pass
    
    f.close()
    
    x_copy,y_copy,p_copy,t_copy = np.array(x_copy),np.array(y_copy),np.array(p_copy),np.array(t_copy)
    x_command,y_command,p_command,t_command = np.array(x_command),np.array(y_command),np.array(p_command),np.array(t_command)

    print len(x_copy),len(y_copy)
    # interpolate to find missing points if necessary
    # assume for right now that we are never missing multiple consecutive points
    sinc = None
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
            sinc = np.concatenate((sinc_neg,np.array([0]),sinc_pos))
            # evaluate convolution at missing point
            missing_x = np.dot(xexpand_copy,sinc)
            missing_y = np.dot(yexpand_copy,sinc)
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
        print len(xintp_copy),len(yintp_copy)

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
    sincxt_copy.stem(np.arange(30),sinc[:30],linestyle='none',marker='o')
    sincyt_copy.stem(np.arange(30),sinc[:30],linestyle='none',marker='o')
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
    xt_copy.set_ylabel('x[n]',fontsize=y_axis_fontsize)
    yt_copy.set_ylabel('y[n]',fontsize=y_axis_fontsize)
    sincxt_copy.set_ylabel('sinc[n]',fontsize=y_axis_fontsize)
    sincyt_copy.set_ylabel('sinc[n]',fontsize=y_axis_fontsize)
    xintpt_copy.set_ylabel('x[n] \n (w/ interpolation)',fontsize=y_axis_fontsize-7)
    yintpt_copy.set_ylabel('y[n] \n (w/ interpolation)',fontsize=y_axis_fontsize-7)

    # save figures
    fig_xt_copy.savefig(path+'figs_raw/'+fname[:len(fname)-4]+'/xintpt_copy_'+fname[:len(fname)-4]+'.png')
    fig_yt_copy.savefig(path+'figs_raw/'+fname[:len(fname)-4]+'/yintpt_copy_'+fname[:len(fname)-4]+'.png')

    # command clocks
    fig_xt_command,fig_yt_command,fig_xy_command = plt.figure(),plt.figure(),plt.figure()
    xt_command,yt_command,xy_command = fig_xt_command.add_subplot(111),fig_yt_command.add_subplot(111),fig_xy_command.add_subplot(111,aspect=1.0)
    xt_command.plot(t_command-t_command[0],x_command)
    yt_command.plot(t_command-t_command[0],y_command)
    xy_command.plot(x_command,y_command)

    # equalize axis scales
    if max(x_command)-min(x_command)>max(y_command)-min(y_command):
        ax_range = max(x_copy)-min(x_command)+20
        xy_command.set_xlim(min(x_command)-10,max(x_command)+10)
        xy_command.set_ylim(min(y_command)-10,min(y_command)+ax_range)
    else:
        ax_range = max(y_command)-min(y_command)+20
        xy_command.set_xlim(min(x_command)-10,min(x_command)+ax_range)
        xy_command.set_ylim(min(y_command)-10,max(y_command)+10)
    plt.axis('equal')

    # set titles and file labels
    #xt_command.set_title('Sample Command Clock',fontsize=title_fontsize)
    #yt_command.set_title('Sample Command Clock',fontsize=title_fontsize)
    #xy_command.set_title('Sample Command Clock',fontsize=title_fontsize)

    fig_xt_command.text(0.99, 0.96, fname[:len(fname)-4],fontsize=10,color='red',va='baseline',ha='right',multialignment='left')
    fig_yt_command.text(0.99, 0.96, fname[:len(fname)-4],fontsize=10,color='red',va='baseline',ha='right',multialignment='left')
    fig_xy_command.text(0.99, 0.96, fname[:len(fname)-4],fontsize=10,color='red',va='baseline',ha='right',multialignment='left')

    if 'YDU' in fname:
        fig_xt_command.text(0.25, 0.955, 'HEALTHY',fontsize=15,color='black',va='baseline',ha='right',multialignment='left')
        fig_yt_command.text(0.25, 0.955, 'HEALTHY',fontsize=15,color='black',va='baseline',ha='right',multialignment='left')
        fig_xy_command.text(0.25, 0.955, 'HEALTHY',fontsize=15,color='black',va='baseline',ha='right',multialignment='left')
    elif 'CIN' in fname:
        fig_xt_command.text(0.25, 0.955, 'IMPAIRED',fontsize=15,color='black',va='baseline',ha='right',multialignment='left')
        fig_yt_command.text(0.25, 0.955, 'IMPAIRED',fontsize=15,color='black',va='baseline',ha='right',multialignment='left')
        fig_xy_command.text(0.25, 0.955, 'IMPAIRED',fontsize=15,color='black',va='baseline',ha='right',multialignment='left')
    else:
        print 'not a valid filename'

    # set axis labels
    xt_command.set_xlabel('time [ms]',fontsize=x_axis_fontsize)
    yt_command.set_xlabel('time [ms]',fontsize=x_axis_fontsize)
    xy_command.set_xlabel('x',fontsize=x_axis_fontsize)

    xt_command.set_ylabel('x',fontsize=y_axis_fontsize)
    yt_command.set_ylabel('y',fontsize=y_axis_fontsize)
    xy_command.set_ylabel('y',fontsize=y_axis_fontsize)

    # save figures
    fig_xt_command.savefig(path+'figs_raw/'+fname[:len(fname)-4]+'/xt_command_'+fname[:len(fname)-4]+'.png')
    fig_yt_command.savefig(path+'figs_raw/'+fname[:len(fname)-4]+'/yt_command_'+fname[:len(fname)-4]+'.png')
    fig_xy_command.savefig(path+'figs_raw/'+fname[:len(fname)-4]+'/xy_command_'+fname[:len(fname)-4]+'.png')

    plt.close('all')
    