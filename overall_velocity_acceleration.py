# plot overall velocity and acceleration of pen

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

    # overall velocity and acceleration
    v_copy,v_command = [],[]
    for v in range(1,len(x_copy)):
        dist = math.sqrt((x_copy[v]-x_copy[v-1])**2+(y_copy[v]-y_copy[v-1])**2)
        time = t_copy[v]-t_copy[v-1]
        v_copy.append(dist/time)
    for d in range(1,len(x_command)):
        dist = math.sqrt((x_command[d]-x_command[d-1])**2+(y_command[d]-y_command[d-1])**2)
        time = t_command[d]-t_command[d-1]
        v_command.append(dist/time)
    a_copy,a_command = [],[]
    for n in range(1,len(v_copy)):
        a_copy.append((v_copy[n]-v_copy[n-1])/(t_copy[n]-t_copy[n-1]))
    for r in range(1,len(v_command)):
        a_command.append((v_command[r]-v_command[r-1])/(t_command[r]-t_command[r-1]))
        
    # copy clocks
    fig_vt_at_copy = plt.figure()
    fig_vt_at_copy.subplots_adjust(right=0.85)
    vt_copy = fig_vt_at_copy.add_subplot(111)
    at_copy = vt_copy.twinx()
    
    #marker_size,marker_edgecolor_v,linecolor_v,linestyle_v,linewidth_v = 5,'darkred','black','k--',2.0
    #vxt_copy.plot(t_copy[1:]-t_copy[0],vx_copy,linestyle_v,lw=linewidth_v,markersize=marker_size)#,markeredgecolor=marker_edgecolor_v)
    
    marker_size,marker_edgecolor_v,linecolor_v,linestyle_v,linewidth_v = 5,'darkred','black','k--',2.0
    vt_copy.plot(t_copy[1:]-t_copy[0],v_copy,linestyle_v,lw=linewidth_v,markersize=marker_size)
    
    marker_edgecolor_a,linecolor_a,linewidth_a = 'darkgreen','green',1.5
    at_copy.plot(t_copy[2:]-t_copy[0],a_copy,linecolor_a,lw=linewidth_a)

    # set axis limits
    vt_copy.set_xlim(right=max(t_copy-t_copy[0]))

    # set titles and file labels
    title_fontsize = 20
    #xt_copy.set_title('Sample Copy Clock',fontsize=title_fontsize)
    #yt_copy.set_title('Sample Copy Clock',fontsize=title_fontsize)
    #xy_copy.set_title('Sample Copy Clock',fontsize=title_fontsize)

    fig_vt_at_copy.text(0.99, 0.96, fname[:len(fname)-4],fontsize=10,color='red',va='baseline',ha='right',multialignment='left')

    if 'YDU' in fname:
        fig_vt_at_copy.text(0.25, 0.955, 'HEALTHY',fontsize=15,color='black',va='baseline',ha='right',multialignment='left')
    elif 'CIN' in fname:
        fig_vt_at_copy.text(0.25, 0.955, 'IMPAIRED',fontsize=15,color='black',va='baseline',ha='right',multialignment='left')
    else:
        print 'not a valid filename'

    # set axis labels
    x_axis_fontsize = 20
    vt_copy.set_xlabel('time [ms]',fontsize=x_axis_fontsize)

    y_axis_fontsize = 20
    vt_copy.set_ylabel(r'$v$',color=linecolor_v,fontsize=y_axis_fontsize)
    at_copy.set_ylabel(r'$a$',color=linecolor_a,fontsize=y_axis_fontsize)
    for v1 in vt_copy.get_yticklabels():
        v1.set_color(linecolor_v)
    for a1 in at_copy.get_yticklabels():
        a1.set_color(linecolor_a)

    # save figures
    fig_vt_at_copy.savefig(path+'figs_raw/'+fname[:len(fname)-4]+'/vt_at_copy_'+fname[:len(fname)-4]+'.png')

    # copy clocks
    fig_vt_at_command = plt.figure()
    fig_vt_at_command.subplots_adjust(right=0.85)
    vt_command = fig_vt_at_command.add_subplot(111)
    at_command = vt_command.twinx()
    
    #marker_size,marker_edgecolor_v,linecolor_v,linestyle_v,linewidth_v = 5,'darkred','black','k--',2.0
    #vxt_command.plot(t_command[1:]-t_command[0],vx_command,linestyle_v,lw=linewidth_v,markersize=marker_size)#,markeredgecolor=marker_edgecolor_v)
    
    marker_size,marker_edgecolor_v,linecolor_v,linestyle_v,linewidth_v = 5,'darkred','black','k--',2.0
    vt_command.plot(t_command[1:]-t_command[0],v_command,linestyle_v,lw=linewidth_v,markersize=marker_size)
    
    marker_edgecolor_a,linecolor_a,linewidth_a = 'darkgreen','green',1.5
    at_command.plot(t_command[2:]-t_command[0],a_command,linecolor_a,lw=linewidth_a)

    # set axis limits
    vt_command.set_xlim(right=max(t_command-t_command[0]))

    # set titles and file labels
    title_fontsize = 20
    #xt_command.set_title('Sample Copy Clock',fontsize=title_fontsize)
    #yt_command.set_title('Sample Copy Clock',fontsize=title_fontsize)
    #xy_command.set_title('Sample Copy Clock',fontsize=title_fontsize)

    fig_vt_at_command.text(0.99, 0.96, fname[:len(fname)-4],fontsize=10,color='red',va='baseline',ha='right',multialignment='left')

    if 'YDU' in fname:
        fig_vt_at_command.text(0.25, 0.955, 'HEALTHY',fontsize=15,color='black',va='baseline',ha='right',multialignment='left')
    elif 'CIN' in fname:
        fig_vt_at_command.text(0.25, 0.955, 'IMPAIRED',fontsize=15,color='black',va='baseline',ha='right',multialignment='left')
    else:
        print 'not a valid filename'

    # set axis labels
    x_axis_fontsize = 20
    vt_command.set_xlabel('time [ms]',fontsize=x_axis_fontsize)

    y_axis_fontsize = 20
    vt_command.set_ylabel(r'$v$',color=linecolor_v,fontsize=y_axis_fontsize)
    at_command.set_ylabel(r'$a$',color=linecolor_a,fontsize=y_axis_fontsize)
    for v1 in vt_command.get_yticklabels():
        v1.set_color(linecolor_v)
    for a1 in at_command.get_yticklabels():
        a1.set_color(linecolor_a)

    # save figures
    fig_vt_at_command.savefig(path+'figs_raw/'+fname[:len(fname)-4]+'/vt_at_command_'+fname[:len(fname)-4]+'.png')
    
    plt.close('all')
