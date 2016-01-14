# -*- coding: utf-8 -*-
import math
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
from pylab import *

# the circles are way off and don’t look helpful at all

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

    # guess true circles
    xcenter_copy,ycenter_copy = np.mean(x_copy),np.mean(y_copy)
    total_dist = 0
    for w in range(len(x_copy)):
        dist_from_center = math.sqrt((x_copy[w]-xcenter_copy)**2+(y_copy[w]-ycenter_copy)**2)
        total_dist += dist_from_center
    avg_dist_from_center = total_dist/len(x_copy)
    radius_copy = avg_dist_from_center
    angles = np.linspace(0,2*math.pi,250)
    xtrue_copy = radius_copy*np.cos(angles)+xcenter_copy
    ytrue_copy = radius_copy*np.sin(angles)+ycenter_copy
    print 'copy clock: radius = ',radius_copy,', center = (',xcenter_copy,',',ycenter_copy,')'
    print xtrue_copy[:5]
    
    xcenter_command,ycenter_command = np.mean(x_command),np.mean(y_command)
    total_dist = 0
    for w in range(len(x_command)):
        dist_from_center = math.sqrt((x_command[w]-xcenter_command)**2+(y_command[w]-ycenter_command)**2)
        total_dist += dist_from_center
    avg_dist_from_center = total_dist/len(x_command)
    radius_command = avg_dist_from_center
    xtrue_command = radius_command*np.cos(angles)+xcenter_command
    ytrue_command = radius_command*np.sin(angles)+ycenter_command
    print 'command clock: radius = ',radius_command,', center = (',xcenter_command,',',ycenter_command,')'
    print ytrue_copy[:5]

    # copy clocks
    fig_xy_copy = plt.figure()
    xy_copy = fig_xy_copy.add_subplot(111,aspect=1.0)
    xy_copy.plot(x_copy,y_copy)
    xy_copy.plot(xtrue_copy,ytrue_copy,'--')

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
    #xy_copy.set_title('Sample Copy Clock',fontsize=title_fontsize)

    fig_xy_copy.text(0.99, 0.96, fname[:len(fname)-4],fontsize=10,color='red',va='baseline',ha='right',multialignment='left')

    if 'YDU' in fname:
        fig_xy_copy.text(0.25, 0.955, 'HEALTHY',fontsize=15,color='black',va='baseline',ha='right',multialignment='left')
    elif 'CIN' in fname:
        fig_xy_copy.text(0.25, 0.955, 'IMPAIRED',fontsize=15,color='black',va='baseline',ha='right',multialignment='left')
    else:
        print 'not a valid filename'

    # set axis labels
    x_axis_fontsize = 20
    xy_copy.set_xlabel('x',fontsize=x_axis_fontsize)

    y_axis_fontsize = 20
    xy_copy.set_ylabel('y',fontsize=y_axis_fontsize)

    # save figures
    fig_xy_copy.savefig(path+'figs_raw/'+fname[:len(fname)-4]+'/xy_copy_'+fname[:len(fname)-4]+'.png')

    # command clocks
    fig_xy_command = plt.figure()
    xy_command = fig_xy_command.add_subplot(111,aspect=1.0)
    xy_command.plot(x_command,y_command)
    xy_command.plot(xtrue_command,ytrue_command,'--')

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
    #xy_command.set_title('Sample Command Clock',fontsize=title_fontsize)

    fig_xy_command.text(0.99, 0.96, fname[:len(fname)-4],fontsize=10,color='red',va='baseline',ha='right',multialignment='left')

    if 'YDU' in fname:
        fig_xy_command.text(0.25, 0.955, 'HEALTHY',fontsize=15,color='black',va='baseline',ha='right',multialignment='left')
    elif 'CIN' in fname:
        fig_xy_command.text(0.25, 0.955, 'IMPAIRED',fontsize=15,color='black',va='baseline',ha='right',multialignment='left')
    else:
        print 'not a valid filename'

    # set axis labels
    xy_command.set_xlabel('x',fontsize=x_axis_fontsize)
    xy_command.set_ylabel('y',fontsize=y_axis_fontsize)

    # save figures
    fig_xy_command.savefig(path+'figs_raw/'+fname[:len(fname)-4]+'/xy_command_'+fname[:len(fname)-4]+'.png')

    plt.close('all')
    