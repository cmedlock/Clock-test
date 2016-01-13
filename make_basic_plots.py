import math
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import numpy.fft
import os
from pylab import *

import clock_test as ct
ct = reload(ct)
    
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
    x_copy,y_copy,p_copy,t_copy,x_command,y_command,p_command,t_command = ct.parse_file(fname,path)
    print '\n',len(x_copy),len(y_copy)
    #print len(x_command),len(y_command),'\n'
    
    # interpolate to find missing points if necessary
    # assume for right now that we are only ever missing 1 point
    
    # for testing purposes, remove a point with an odd index from x_copy
    #x_copy[5] = -5 # a point on the left
    x_copy[113] = -5 # a point in the middle
    x_copy[200] = -5
    #x_copy[93] = -5
    #x_copy[225] = -5 # a point on the right
    
    x_copy = ct.interpolate(x_copy)
    print x_copy[5]
    #x_command,y_command = ct.linear_interpolate(x_command,y_command)

    x_copy,y_copy,p_copy,t_copy = np.array(x_copy),np.array(y_copy),np.array(p_copy),np.array(t_copy)
    #x_command,y_command,p_command,t_command = np.array(x_command),np.array(y_command),np.array(p_command),np.array(t_command)
    
    # draw clocks
    ct.draw_clock(x_copy,y_copy,'xy_copy',fname,path)
    #ct.draw_clock(x_command,y_command,'xy_command',fname,path)

    # plot x[t], y[t], v[t], and a[t]
    vx_copy,vy_copy,ax_copy,ay_copy = ct.deriv_doublederiv(x_copy,y_copy,t_copy)
    #vx_command,vy_command,ax_command,ay_command = ct.deriv_doublederiv(x_command,y_command,t_command)
    # copy clocks
    ct.plot_xyt_other(x_copy,t_copy-t_copy[0],'x','time [ms]',vx_copy,t_copy[1:]-t_copy[0],'v_x','time [ms]',False,'vxt_copy',fname,path)
    """
    ct.plot_xyt_other(x_copy,t_copy-t_copy[0],'x','time [ms]',ax_copy,t_copy[2:]-t_copy[0],'a_x','time [ms]',False,'axt_copy',fname,path)
    ct.plot_xyt_other(y_copy,t_copy-t_copy[0],'y','time [ms]',vy_copy,t_copy[1:]-t_copy[0],'v_y','time [ms]',False,'vyt_copy',fname,path)
    ct.plot_xyt_other(y_copy,t_copy-t_copy[0],'y','time [ms]',ay_copy,t_copy[2:]-t_copy[0],'a_y','time [ms]',False,'ayt_copy',fname,path)
    # command clocks
    ct.plot_xyt_other(x_command,t_command-t_command[0],'x','time [ms]',vx_command,t_command[1:]-t_command[0],'v_x','time [ms]',False,'vxt_command',fname,path)
    ct.plot_xyt_other(x_command,t_command-t_command[0],'x','time [ms]',ax_command,t_command[2:]-t_command[0],'a_x','time [ms]',False,'axt_command',fname,path)
    ct.plot_xyt_other(y_command,t_command-t_command[0],'y','time [ms]',vy_command,t_command[1:]-t_command[0],'v_y','time [ms]',False,'vyt_command',fname,path)
    ct.plot_xyt_other(y_command,t_command-t_command[0],'y','time [ms]',ay_command,t_command[2:]-t_command[0],'a_y','time [ms]',False,'ayt_command',fname,path)

    # plot spectrum of x[t] and y[t]
    
    # plot pressure
    """
print 'Done'
