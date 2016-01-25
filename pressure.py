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
            x_copy.append(xcoord)
            y_copy.append(ycoord)
            p_copy.append(pressure)
            t_copy.append(timestamp)
        elif clock_type=='COMMAND':
            x_command.append(xcoord)
            y_command.append(ycoord)
            p_command.append(pressure)
            t_command.append(timestamp)
        else:
            print 'not a valid clock type'
    
    f.close()

    x_copy,y_copy,p_copy,t_copy = np.array(x_copy),np.array(y_copy),np.array(p_copy),np.array(t_copy)
    x_command,y_command,p_command,t_command = np.array(x_command),np.array(y_command),np.array(p_command),np.array(t_command)

    # copy clocks
    fig_pt_copy = plt.figure()
    fig_pt_copy.subplots_adjust(hspace=0.3)
    pt_copy = fig_pt_copy.add_subplot(111)
    linecolor_p = 'red'
    pt_copy.plot(t_copy-t_copy[0],p_copy,linecolor_p)
    
    # set axis limits
    pt_copy.set_xlim(right=max(t_copy-t_copy[0]))
    
    # set titles and file labels
    title_fontsize = 20
    #xt_copy.set_title('Sample Copy Clock',fontsize=title_fontsize)
    #yt_copy.set_title('Sample Copy Clock',fontsize=title_fontsize)
    #xy_copy.set_title('Sample Copy Clock',fontsize=title_fontsize)

    fig_pt_copy.text(0.99, 0.96,fname[:len(fname)-4],fontsize=10,color='red',va='baseline',ha='right',multialignment='left')

    if 'YDU' in fname:
        fig_pt_copy.text(0.25, 0.955, 'HEALTHY',fontsize=15,color='black',va='baseline',ha='right',multialignment='left')
    elif 'CIN' in fname:
        fig_pt_copy.text(0.25, 0.955, 'IMPAIRED',fontsize=15,color='black',va='baseline',ha='right',multialignment='left')
    else:
        print 'not a valid filename'
        
    # set axis labels
    x_axis_fontsize = 20
    pt_copy.set_xlabel('time [ms]',fontsize=x_axis_fontsize)

    y_axis_fontsize = 25
    pt_copy.set_ylabel('pressure',color=linecolor_p,fontsize=y_axis_fontsize)
    for p1 in pt_copy.get_yticklabels():
        p1.set_color(linecolor_p)
        
    # save figures
    fig_pt_copy.savefig(path+'figs_raw/'+fname[:len(fname)-4]+'/pt_copy_'+fname[:len(fname)-4]+'.png')

    # command clocks
    fig_pt_command = plt.figure()
    fig_pt_command.subplots_adjust(hspace=0.3)
    pt_command = fig_pt_command.add_subplot(111)
    pt_command.plot(t_command-t_command[0],p_command,linecolor_p)

    # set axis limits
    pt_command.set_xlim(right=max(t_command-t_command[0]))
    
    # set titles and file labels
    #xt_command.set_title('Sample Command Clock',fontsize=title_fontsize)
    #yt_command.set_title('Sample Command Clock',fontsize=title_fontsize)
    #xy_command.set_title('Sample Command Clock',fontsize=title_fontsize)

    fig_pt_command.text(0.99, 0.96, fname[:len(fname)-4],fontsize=10,color='red',va='baseline',ha='right',multialignment='left')

    if 'YDU' in fname:
        fig_pt_command.text(0.25, 0.955, 'HEALTHY',fontsize=15,color='black',va='baseline',ha='right',multialignment='left')
    elif 'CIN' in fname:
        fig_pt_command.text(0.25, 0.955, 'IMPAIRED',fontsize=15,color='black',va='baseline',ha='right',multialignment='left')
    else:
        print 'not a valid filename'
        
    # set axis labels
    pt_command.set_xlabel('time [ms]',fontsize=x_axis_fontsize)

    pt_command.set_ylabel('pressure',color=linecolor_p,fontsize=y_axis_fontsize)
    for p1 in pt_command.get_yticklabels():
        p1.set_color(linecolor_p)

    # save figures
    fig_pt_command.savefig(path+'figs_raw/'+fname[:len(fname)-4]+'/pt_command_'+fname[:len(fname)-4]+'.png')

    plt.close('all')
    