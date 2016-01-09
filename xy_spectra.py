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

for fname in dirs:
    if 'Scored' not in fname:
        continue
    print 'reading file ',fname,'...'
    f = open(path+fname)
    data = f.readlines()

    if not os.path.exists(path+'figs_raw/'+fname[:len(fname)-4]):
        os.makedirs(path+'figs_raw/'+fname[:len(fname)-4])

    # copy or command clock?
    clock_type = ''
    x_copy,y_copy,t_copy = [],[],[]
    x_command,y_command,t_command = [],[],[]

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
        timestamp = double(line[7])
        if clock_type=='COPY':
            x_copy.append(xcoord)
            y_copy.append(ycoord)
            t_copy.append(timestamp)
        elif clock_type=='COMMAND':
            x_command.append(xcoord)
            y_command.append(ycoord)
            t_command.append(timestamp)
        else:
            print 'not a valid clock type'
    
    f.close()

    x_copy,y_copy,t_copy = np.array(x_copy),np.array(y_copy),np.array(t_copy)
    x_command,y_command,t_command = np.array(x_command),np.array(y_command),np.array(t_command)

    # dft's
    dft_size = 750
    dftx_copy,dfty_copy = np.fft.fft(x_copy,n=dft_size),np.fft.fft(y_copy,n=dft_size)
    dftx_command,dfty_command = np.fft.fft(x_command,n=dft_size),np.fft.fft(y_command,n=dft_size)
    k = np.arange(dft_size)
    freq = 2*np.pi/dft_size*k
    #freq = np.fft.fftfreq(n=1024,d=1/(2*np.pi)) # NB: centered around 0
    
    # copy clocks
    fig_xt_copy,fig_yt_copy = plt.figure(),plt.figure()
    fig_xt_copy.subplots_adjust(hspace=0.3)
    fig_yt_copy.subplots_adjust(hspace=0.3)
    xt_copy,yt_copy = fig_xt_copy.add_subplot(211),fig_yt_copy.add_subplot(211)
    xt_copy.plot(t_copy-t_copy[0],x_copy)
    yt_copy.plot(t_copy-t_copy[0],y_copy)
    
    dftxk_copy,dftyk_copy = fig_xt_copy.add_subplot(212),fig_yt_copy.add_subplot(212)    
    linecolor_dft = 'red'
    dftxk_copy.plot(k,np.abs(dftx_copy),linecolor_dft)
    dftyk_copy.plot(k,np.abs(dfty_copy),linecolor_dft)
    dftxk_copy.set_yscale('log')
    dftyk_copy.set_yscale('log')

    # set titles and file labels
    title_fontsize = 20
    #xt_copy.set_title('Sample Copy Clock',fontsize=title_fontsize)
    #yt_copy.set_title('Sample Copy Clock',fontsize=title_fontsize)
    #xy_copy.set_title('Sample Copy Clock',fontsize=title_fontsize)

    fig_xt_copy.text(0.99, 0.96,fname[:len(fname)-4],fontsize=10,color='red',va='baseline',ha='right',multialignment='left')
    fig_yt_copy.text(0.99, 0.96,fname[:len(fname)-4],fontsize=10,color='red',va='baseline',ha='right',multialignment='left')

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
    xt_copy.set_xlabel('time [ms]',fontsize=x_axis_fontsize)
    yt_copy.set_xlabel('time [ms]',fontsize=x_axis_fontsize)
    dftxk_copy.set_xlabel(r'$k$',fontsize=x_axis_fontsize)
    dftyk_copy.set_xlabel(r'$k$',fontsize=x_axis_fontsize)

    y_axis_fontsize = 25
    linecolor_xy = 'blue'
    xt_copy.set_ylabel(r'$x$',color=linecolor_xy,fontsize=y_axis_fontsize)
    yt_copy.set_ylabel(r'$y$',color=linecolor_xy,fontsize=y_axis_fontsize)
    for x1 in xt_copy.get_yticklabels():
        x1.set_color(linecolor_xy)
    for y1 in yt_copy.get_yticklabels():
        y1.set_color(linecolor_xy)
        
    dftxk_copy.set_ylabel(r'$|X[k]|$',color=linecolor_dft,fontsize=y_axis_fontsize)
    dftyk_copy.set_ylabel(r'$|Y[k]|$',color=linecolor_dft,fontsize=y_axis_fontsize)
    for v1 in dftxk_copy.get_yticklabels():
        v1.set_color(linecolor_dft)
    for v2 in dftyk_copy.get_yticklabels():
        v2.set_color(linecolor_dft)
    
    # save figures
    fig_xt_copy.savefig(path+'figs_raw/'+fname[:len(fname)-4]+'/dftx_copy_'+fname[:len(fname)-4]+'.png')
    fig_yt_copy.savefig(path+'figs_raw/'+fname[:len(fname)-4]+'/dfty_copy_'+fname[:len(fname)-4]+'.png')

    # command clocks
    fig_xt_command,fig_yt_command = plt.figure(),plt.figure()
    fig_xt_command.subplots_adjust(hspace=0.3)
    fig_yt_command.subplots_adjust(hspace=0.3)
    xt_command,yt_command = fig_xt_command.add_subplot(211),fig_yt_command.add_subplot(211)
    xt_command.plot(t_command-t_command[0],x_command)
    yt_command.plot(t_command-t_command[0],y_command)

    dftxk_command,dftyk_command = fig_xt_command.add_subplot(212),fig_yt_command.add_subplot(212)    
    linecolor_dft = 'red'
    dftxk_command.plot(k,np.abs(dftx_command),linecolor_dft)
    dftyk_command.plot(k,np.abs(dfty_command),linecolor_dft)
    dftxk_command.set_yscale('log')
    dftyk_command.set_yscale('log')

    # set titles and file labels
    #xt_command.set_title('Sample Command Clock',fontsize=title_fontsize)
    #yt_command.set_title('Sample Command Clock',fontsize=title_fontsize)
    #xy_command.set_title('Sample Command Clock',fontsize=title_fontsize)

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
    xt_command.set_xlabel('time [ms]',fontsize=x_axis_fontsize)
    yt_command.set_xlabel('time [ms]',fontsize=x_axis_fontsize)
    dftxk_command.set_xlabel(r'$k$',fontsize=x_axis_fontsize)
    dftyk_command.set_xlabel(r'$k$',fontsize=x_axis_fontsize)

    xt_command.set_ylabel(r'$x$',color=linecolor_xy,fontsize=y_axis_fontsize)
    yt_command.set_ylabel(r'$y$',color=linecolor_xy,fontsize=y_axis_fontsize)
    for x1 in xt_command.get_yticklabels():
        x1.set_color(linecolor_xy)
    for y1 in yt_command.get_yticklabels():
        y1.set_color(linecolor_xy)

    dftxk_command.set_ylabel(r'$|X[k]|$',color=linecolor_dft,fontsize=y_axis_fontsize)
    dftyk_command.set_ylabel(r'$|Y[k]|$',color=linecolor_dft,fontsize=y_axis_fontsize)
    for v1 in dftxk_command.get_yticklabels():
        v1.set_color(linecolor_dft)
    for v2 in dftyk_command.get_yticklabels():
        v2.set_color(linecolor_dft)

    # save figures
    fig_xt_command.savefig(path+'figs_raw/'+fname[:len(fname)-4]+'/dftx_command_'+fname[:len(fname)-4]+'.png')
    fig_yt_command.savefig(path+'figs_raw/'+fname[:len(fname)-4]+'/dfty_command_'+fname[:len(fname)-4]+'.png')

    plt.close('all')
    