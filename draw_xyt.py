import math
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
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

    # copy clocks
    fig_xyt_copy,fig_xy_copy = plt.figure(),plt.figure()
    xt_copy,yt_copy,xy_copy = fig_xyt_copy.add_subplot(211),fig_xyt_copy.add_subplot(212),fig_xy_copy.add_subplot(111,aspect=1.0)
    xt_copy.plot(t_copy-t_copy[0],x_copy)
    yt_copy.plot(t_copy-t_copy[0],y_copy)
    xy_copy.plot(x_copy,y_copy)

    # equalize axis scales
    if max(x_copy)-min(x_copy)>max(y_copy)-min(y_copy):
        ax_range = max(x_copy)-min(x_copy)+20
        xy_copy.set_xlim(min(x_copy)-10,max(x_copy)+10)
        xy_copy.set_ylim(min(y_copy)-10,min(y_copy)+ax_range)
    else:
        ax_range = max(y_copy)-min(y_copy)+20
        xy_copy.set_xlim(min(x_copy)-10,min(x_copy)+ax_range)
        xy_copy.set_ylim(min(y_copy)-10,max(y_copy)+10)
    plt.axis('equal')
    
    # set titles and file labels
    title_fontsize = 20
    #xt_copy.set_title('Sample Copy Clock',fontsize=title_fontsize)
    #yt_copy.set_title('Sample Copy Clock',fontsize=title_fontsize)
    #xy_copy.set_title('Sample Copy Clock',fontsize=title_fontsize)

    fig_xyt_copy.text(0.99, 0.96, fname[:len(fname)-4],fontsize=10,color='red',va='baseline',ha='right',multialignment='left')
    fig_xy_copy.text(0.99, 0.96, fname[:len(fname)-4],fontsize=10,color='red',va='baseline',ha='right',multialignment='left')

    if 'YDU' in fname:
        fig_xyt_copy.text(0.25, 0.955, 'HEALTHY',fontsize=15,color='black',va='baseline',ha='right',multialignment='left')
        fig_xy_copy.text(0.25, 0.955, 'HEALTHY',fontsize=15,color='black',va='baseline',ha='right',multialignment='left')
    elif 'CIN' in fname:
        fig_xyt_copy.text(0.25, 0.955, 'IMPAIRED',fontsize=15,color='black',va='baseline',ha='right',multialignment='left')
        fig_xy_copy.text(0.25, 0.955, 'IMPAIRED',fontsize=15,color='black',va='baseline',ha='right',multialignment='left')
    else:
        print 'not a valid filename'

    # set axis labels
    x_axis_fontsize = 20
    yt_copy.set_xlabel('time [ms]',fontsize=x_axis_fontsize)
    xy_copy.set_xlabel('x',fontsize=x_axis_fontsize)

    y_axis_fontsize = 20
    xt_copy.set_ylabel('x',fontsize=y_axis_fontsize)
    yt_copy.set_ylabel('y',fontsize=y_axis_fontsize)
    xy_copy.set_ylabel('y',fontsize=y_axis_fontsize)

    # save figures
    fig_xyt_copy.savefig(path+'figs_raw/'+fname[:len(fname)-4]+'/xyt_copy_'+fname[:len(fname)-4]+'.png')
    fig_xy_copy.savefig(path+'figs_raw/'+fname[:len(fname)-4]+'/xy_copy_'+fname[:len(fname)-4]+'.png')

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
    