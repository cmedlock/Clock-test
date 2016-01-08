import math
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
from pylab import *

path = '/Users/cmedlock/Documents/DSP UROP/all_data/'
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

    # velocities and accelerations
    vx_copy,vy_copy = [],[]
    for d in range(1,len(x_copy)):
        vx = (x_copy[d]-x_copy[d-1])/(t_copy[d]-t_copy[d-1])
        vy = (y_copy[d]-y_copy[d-1])/(t_copy[d]-t_copy[d-1])
        vx_copy.append(vx)
        vy_copy.append(vy)
    vx_command,vy_command = [],[]
    for a in range(1,len(x_command)):
        vx = (x_command[a]-x_command[a-1])/(t_command[a]-t_command[a-1])
        vy = (y_command[a]-y_command[a-1])/(t_command[a]-t_command[a-1])
        vx_command.append(vx)
        vy_command.append(vy)
    ax_copy,ay_copy = [],[]
    for b in range(1,len(vx_copy)):
        ax = (vx_copy[b]-vx_copy[b-1])/(t_copy[b]-t_copy[b-1])
        ay = (vy_copy[b]-vy_copy[b-1])/(t_copy[b]-t_copy[b-1])
        ax_copy.append(ax)
        ay_copy.append(ay)
    ax_command,ay_command = [],[]
    for c in range(1,len(vx_command)):
        ax = (vx_command[c]-vx_command[c-1])/(t_command[c]-t_command[c-1])
        ay = (vy_command[c]-vy_command[c-1])/(t_command[c]-t_command[c-1])
        ax_command.append(ax)
        ay_command.append(ay)
        
    # copy clocks
    fig_vxt_copy,fig_vyt_copy = plt.figure(),plt.figure()
    fig_vxt_copy.subplots_adjust(right=0.85)
    fig_vyt_copy.subplots_adjust(right=0.85)
    vxt_copy,vyt_copy = fig_vxt_copy.add_subplot(111),fig_vyt_copy.add_subplot(111)    
    marker_size,marker_edgecolor_v,linecolor_v = 5,'darkred','red'
    vxt_copy.plot(t_copy[1:]-t_copy[0],vx_copy,linecolor_v,markersize=marker_size)#,markeredgecolor=marker_edgecolor_v)
    vyt_copy.plot(t_copy[1:]-t_copy[0],vy_copy,linecolor_v,markersize=marker_size)#,markeredgecolor=marker_edgecolor_v)

    axt_copy,ayt_copy = vxt_copy.twinx(),vyt_copy.twinx()
    marker_edgecolor_a,linecolor_a = 'darkgreen','green'
    axt_copy.plot(t_copy[2:]-t_copy[0],ax_copy,linecolor_a,markersize=marker_size)#,markeredgecolor=marker_edgecolor_a)
    ayt_copy.plot(t_copy[2:]-t_copy[0],ay_copy,linecolor_a,markersize=marker_size)#,markeredgecolor=marker_edgecolor_a)
    
    title_fontsize = 20
    #xt_copy.set_title('Sample Copy Clock',fontsize=title_fontsize)
    #yt_copy.set_title('Sample Copy Clock',fontsize=title_fontsize)
    #xy_copy.set_title('Sample Copy Clock',fontsize=title_fontsize)

    fig_vxt_copy.text(0.99, 0.96, fname[:len(fname)-4],fontsize=10,color='red',va='baseline',ha='right',multialignment='left')
    fig_vyt_copy.text(0.99, 0.96, fname[:len(fname)-4],fontsize=10,color='red',va='baseline',ha='right',multialignment='left')

    x_axis_fontsize = 20
    vxt_copy.set_xlabel('time [ms]',fontsize=x_axis_fontsize)
    vyt_copy.set_xlabel('time [ms]',fontsize=x_axis_fontsize)

    y_axis_fontsize = 25
    vxt_copy.set_ylabel(r'$v_x$',color=linecolor_v,fontsize=y_axis_fontsize)
    vyt_copy.set_ylabel(r'$v_y$',color=linecolor_v,fontsize=y_axis_fontsize)
    axt_copy.set_ylabel(r'$a_x$',color=linecolor_a,fontsize=y_axis_fontsize)
    ayt_copy.set_ylabel(r'$a_y$',color=linecolor_a,fontsize=y_axis_fontsize)
    for v1 in vxt_copy.get_yticklabels():
        v1.set_color(linecolor_v)
    for v2 in vyt_copy.get_yticklabels():
        v2.set_color(linecolor_v)
    for m1 in axt_copy.get_yticklabels():
        m1.set_color(linecolor_a)
    for m2 in ayt_copy.get_yticklabels():
        m2.set_color(linecolor_a)
    
    fig_vxt_copy.savefig(path+'figs_raw/'+fname[:len(fname)-4]+'/vxt_axt_copy_'+fname[:len(fname)-4]+'.png')
    fig_vyt_copy.savefig(path+'figs_raw/'+fname[:len(fname)-4]+'/vyt_ayt_copy_'+fname[:len(fname)-4]+'.png')

    # command clocks
    fig_vxt_command,fig_vyt_command = plt.figure(),plt.figure()
    fig_vxt_command.subplots_adjust(right=0.85)
    fig_vyt_command.subplots_adjust(right=0.85)
    vxt_command,vyt_command = fig_vxt_command.add_subplot(111),fig_vyt_command.add_subplot(111)
    vxt_command.plot(t_command[1:]-t_command[0],vx_command,linecolor_v,markersize=marker_size)#,markeredgecolor=marker_edgecolor_v)
    vyt_command.plot(t_command[1:]-t_command[0],vy_command,linecolor_v,markersize=marker_size)#,markeredgecolor=marker_edgecolor_v)

    axt_command,ayt_command = vxt_command.twinx(),vyt_command.twinx()
    axt_command.plot(t_command[2:]-t_command[0],ax_command,linecolor_a,markersize=marker_size)#,markeredgecolor=marker_edgecolor_a)
    ayt_command.plot(t_command[2:]-t_command[0],ay_command,linecolor_a,markersize=marker_size)#,markeredgecolor=marker_edgecolor_a)

    #xt_command.set_title('Sample Command Clock',fontsize=title_fontsize)
    #yt_command.set_title('Sample Command Clock',fontsize=title_fontsize)
    #xy_command.set_title('Sample Command Clock',fontsize=title_fontsize)

    fig_vxt_command.text(0.99, 0.96, fname[:len(fname)-4],fontsize=10,color='red',va='baseline',ha='right',multialignment='left')
    fig_vyt_command.text(0.99, 0.96, fname[:len(fname)-4],fontsize=10,color='red',va='baseline',ha='right',multialignment='left')

    vxt_command.set_xlabel('time [ms]',fontsize=x_axis_fontsize)
    vyt_command.set_xlabel('time [ms]',fontsize=x_axis_fontsize)

    vxt_command.set_ylabel(r'$v_x$',color=linecolor_v,fontsize=y_axis_fontsize)
    vyt_command.set_ylabel(r'$v_y$',color=linecolor_v,fontsize=y_axis_fontsize)
    axt_command.set_ylabel(r'$a_x$',color=linecolor_a,fontsize=y_axis_fontsize)
    ayt_command.set_ylabel(r'$a_y$',color=linecolor_a,fontsize=y_axis_fontsize)
    for v1 in vxt_command.get_yticklabels():
        v1.set_color(linecolor_v)
    for v2 in vyt_command.get_yticklabels():
        v2.set_color(linecolor_v)
    for m1 in axt_command.get_yticklabels():
        m1.set_color(linecolor_a)
    for m2 in ayt_command.get_yticklabels():
        m2.set_color(linecolor_a)

    fig_vxt_command.savefig(path+'figs_raw/'+fname[:len(fname)-4]+'/vxt_axt_command_'+fname[:len(fname)-4]+'.png')
    fig_vyt_command.savefig(path+'figs_raw/'+fname[:len(fname)-4]+'/vyt_ayt_command_'+fname[:len(fname)-4]+'.png')

    plt.close('all')
    