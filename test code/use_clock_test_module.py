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
    
    # interpolate to find missing points if necessary
    x_copy,y_copy,p_copy = ct.interpolate(x_copy),ct.interpolate(y_copy),ct.interpolate(p_copy)
    x_command,y_command,p_command = ct.interpolate(x_command),ct.interpolate(y_command),ct.interpolate(p_command)

    # convert to np.arrays before plotting
    x_copy,y_copy,p_copy,t_copy = np.array(x_copy),np.array(y_copy),np.array(p_copy),np.array(t_copy)
    x_command,y_command,p_command,t_command = np.array(x_command),np.array(y_command),np.array(p_command),np.array(t_command)
    
    # draw clocks
    ct.draw_clock(x_copy,y_copy,'x','y','xy_copy',fname,path)
    ct.draw_clock(x_command,y_command,'x','y','xy_command',fname,path)

    # plot x[t], y[t], v[t], and a[t]
    vx_copy,vy_copy,ax_copy,ay_copy = ct.deriv_doublederiv(x_copy,y_copy,t_copy)
    vx_command,vy_command,ax_command,ay_command = ct.deriv_doublederiv(x_command,y_command,t_command)
    # copy clocks
    #ct.plot_xyt_other(x_copy,t_copy-t_copy[0],'x','time [ms]',vx_copy,t_copy[1:]-t_copy[0],'v_x','time [ms]',False,'vxt_copy',fname,path)
    #ct.plot_xyt_other(x_copy,t_copy-t_copy[0],'x','time [ms]',ax_copy,t_copy[2:]-t_copy[0],'a_x','time [ms]',False,'axt_copy',fname,path)
    #ct.plot_xyt_other(y_copy,t_copy-t_copy[0],'y','time [ms]',vy_copy,t_copy[1:]-t_copy[0],'v_y','time [ms]',False,'vyt_copy',fname,path)
    #ct.plot_xyt_other(y_copy,t_copy-t_copy[0],'y','time [ms]',ay_copy,t_copy[2:]-t_copy[0],'a_y','time [ms]',False,'ayt_copy',fname,path)
    # command clocks
    #ct.plot_xyt_other(x_command,t_command-t_command[0],'x','time [ms]',vx_command,t_command[1:]-t_command[0],'v_x','time [ms]',False,'vxt_command',fname,path)
    #ct.plot_xyt_other(x_command,t_command-t_command[0],'x','time [ms]',ax_command,t_command[2:]-t_command[0],'a_x','time [ms]',False,'axt_command',fname,path)
    #ct.plot_xyt_other(y_command,t_command-t_command[0],'y','time [ms]',vy_command,t_command[1:]-t_command[0],'v_y','time [ms]',False,'vyt_command',fname,path)
    #ct.plot_xyt_other(y_command,t_command-t_command[0],'y','time [ms]',ay_command,t_command[2:]-t_command[0],'a_y','time [ms]',False,'ayt_command',fname,path)

    # plot spectrum of x[t] and y[t]
    # subtract mean values (zm = zero mean)
    x_copy_zm,y_copy_zm = x_copy-mean(x_copy),y_copy-mean(y_copy)
    x_command_zm,y_command_zm = x_command-mean(x_command),y_command-mean(y_command)    
    # dft's
    dft_size = 750
    dftx_copy,dfty_copy = np.fft.fft(x_copy_zm,n=dft_size),np.fft.fft(y_copy_zm,n=dft_size)
    dftx_command,dfty_command = np.fft.fft(x_command_zm,n=dft_size),np.fft.fft(y_command_zm,n=dft_size)
    k = np.arange(dft_size)
    freq = 2*np.pi/dft_size*k
    # copy clocks
    ct.plot_xyt_other(x_copy,t_copy-t_copy[0],'x','time [ms]',np.abs(dftx_copy)[:dft_size/2],k[:dft_size/2],'|X[k]|','k',True,'dftx_copy',fname,path)
    ct.plot_xyt_other(y_copy,t_copy-t_copy[0],'y','time [ms]',np.abs(dfty_copy)[:dft_size/2],k[:dft_size/2],'|Y[k]|','k',True,'dfty_copy',fname,path)
    
    # plot pressure
    #ct.draw_clock(t_copy-t_copy[0],p_copy,'time [ms]','pressure','pressure_copy',fname,path)
    #ct.draw_clock(t_command-t_command[0],p_command,'time [ms]','pressure','pressure_command',fname,path)
    
print 'Done'
