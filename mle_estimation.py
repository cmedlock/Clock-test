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

# quantities to compare between healthy and impaired patients
mean_velx_copy,mean_accx_copy,stddev_velx_copy,stddev_accx_copy = [],[],[],[]
mean_vely_copy,mean_accy_copy,stddev_vely_copy,stddev_accy_copy = [],[],[],[]
mean_pressure_copy,stddev_pressure_copy = [],[]
relsidelobe_x_copy,relsidelobe_y_copy = [],[]
mean_velx_command,mean_accx_command,stddev_velx_command,stddev_accx_command = [],[],[],[]
mean_vely_command,mean_accy_command,stddev_vely_command,stddev_accy_command = [],[],[],[]
mean_pressure_command,stddev_pressure_command = [],[]
relsidelobe_x_command,relsidelobe_y_command = [],[]

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

    # get vx[t],vy[t] and ax[t],ay[t]
    vx_copy,vy_copy,ax_copy,ay_copy = ct.deriv_doublederiv(x_copy,y_copy,t_copy)
    vx_command,vy_command,ax_command,ay_command = ct.deriv_doublederiv(x_command,y_command,t_command)
    # record mean velocities and accelerations for comparison between healthy and impaired patients
    # copy clocks
    mean_velx_copy.append(np.mean(vx_copy))
    mean_vely_copy.append(np.mean(vy_copy))
    mean_accx_copy.append(np.mean(ax_copy))
    mean_accy_copy.append(np.mean(ay_copy))
    stddev_velx_copy.append(np.std(vx_copy))
    stddev_vely_copy.append(np.std(vy_copy))
    stddev_accx_copy.append(np.std(ax_copy))
    stddev_accy_copy.append(np.std(ay_copy))
    # command clocks
    mean_velx_command.append(np.mean(vx_command))
    mean_vely_command.append(np.mean(vy_command))
    mean_accx_command.append(np.mean(ax_command))
    mean_accy_command.append(np.mean(ay_command))
    stddev_velx_command.append(np.std(vx_command))
    stddev_vely_command.append(np.std(vy_command))
    stddev_accx_command.append(np.std(ax_command))
    stddev_accy_command.append(np.std(ay_command))

    # get spectra of x[t] and y[t]
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
    abs_dftx_copy = np.abs(dftx_copy)[:dft_size/2] # NB: these are only the positive frequencies
    relsidelobe_x_copy.append(ct.find_relsidelobe_amplitude(abs_dftx_copy))
