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

if not os.path.exists(path+'compare_healthy_impaired'):
    os.makedirs(path+'compare_healthy_impaired')

# quantities to compare the drawings of between healthy and impaired patients
mean_velx_copy,mean_accx_copy,stddev_velx_copy,stddev_accx_copy = [],[],[],[]
mean_vely_copy,mean_accy_copy,stddev_vely_copy,stddev_accy_copy = [],[],[],[]
mean_pressure_copy,stddev_pressure_copy = [],[]
relsidelobe_x_copy,relsidelobe_y_copy = [],[]
mean_velx_command,mean_accx_command,stddev_velx_command,stddev_accx_command = [],[],[],[]
mean_vely_command,mean_accy_command,stddev_vely_command,stddev_accy_command = [],[],[],[]
mean_pressure_command,stddev_pressure_command = [],[]
relsidelobe_x_command,relsidelobe_y_command = [],[]

for fname in dirs:
    
    if 'Scored' not in fname:
        continue
    ftype = ''

    if 'YDU' in fname:
        ftype = 'healthy'
    elif 'CIN' in fname:
        ftype = 'impaired'
    else:
        print 'not a valid file name'
        
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
    mean_velx_copy.append((ftype,np.mean(vx_copy)))
    mean_vely_copy.append((ftype,np.mean(vy_copy)))
    mean_accx_copy.append((ftype,np.mean(ax_copy)))
    mean_accy_copy.append((ftype,np.mean(ay_copy)))
    stddev_velx_copy.append((ftype,np.std(vx_copy)))
    stddev_vely_copy.append((ftype,np.std(vy_copy)))
    stddev_accx_copy.append((ftype,np.std(ax_copy)))
    stddev_accy_copy.append((ftype,np.std(ay_copy)))
    # command clocks
    mean_velx_command.append((ftype,np.mean(vx_command)))
    mean_vely_command.append((ftype,np.mean(vy_command)))
    mean_accx_command.append((ftype,np.mean(ax_command)))
    mean_accy_command.append((ftype,np.mean(ay_command)))
    stddev_velx_command.append((ftype,np.std(vx_command)))
    stddev_vely_command.append((ftype,np.std(vy_command)))
    stddev_accx_command.append((ftype,np.std(ax_command)))
    stddev_accy_command.append((ftype,np.std(ay_command)))

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
    relsidelobe_x_copy.append((ftype,ct.find_relsidelobe_amplitude(dftx_copy)))
    relsidelobe_y_copy.append((ftype,ct.find_relsidelobe_amplitude(dfty_copy)))
    # command clocks
    relsidelobe_x_command.append((ftype,ct.find_relsidelobe_amplitude(dftx_command)))
    relsidelobe_y_command.append((ftype,ct.find_relsidelobe_amplitude(dfty_command)))
    
    # pressure
    mean_pressure_copy.append((ftype,np.mean(p_copy)))
    mean_pressure_command.append((ftype,np.mean(p_command)))
    stddev_pressure_copy.append((ftype,np.std(p_copy)))
    stddev_pressure_command.append((ftype,np.std(p_command)))

# compare quantities from drawings of healthy and impaired patients
# copy clocks
binedges = ct.get_bins(mean_velx_copy,nbins=10)
ct.make_hist([elt[1] for elt in mean_velx_copy if elt[0]=='healthy'],
             [elt[1] for elt in mean_velx_copy if elt[0]=='impaired'],
             binedges,'Mean x Velocity','mean_velx_copy',path)
binedges = ct.get_bins(mean_vely_copy,nbins=10)
ct.make_hist([elt[1] for elt in mean_vely_copy if elt[0]=='healthy'],
             [elt[1] for elt in mean_vely_copy if elt[0]=='impaired'],
             binedges,'Mean y Velocity','mean_vely_copy',path)
binedges = ct.get_bins(stddev_velx_copy,nbins=10)
ct.make_hist([elt[1] for elt in stddev_velx_copy if elt[0]=='healthy'],
             [elt[1] for elt in stddev_velx_copy if elt[0]=='impaired'],
             binedges,'Standard Deviation of x Velocity','stddev_velx_copy',path)
binedges = ct.get_bins(stddev_vely_copy,nbins=10)
ct.make_hist([elt[1] for elt in stddev_vely_copy if elt[0]=='healthy'],
             [elt[1] for elt in stddev_vely_copy if elt[0]=='impaired'],
             binedges,'Standard Deviation of y Velocity','stddev_vely_copy',path)
binedges = ct.get_bins(mean_accx_copy,nbins=10)
ct.make_hist([elt[1] for elt in mean_accx_copy if elt[0]=='healthy'],
             [elt[1] for elt in mean_accx_copy if elt[0]=='impaired'],
             binedges,'Mean x Acceleration','mean_accx_copy',path)
binedges = ct.get_bins(mean_accy_copy,nbins=10)
ct.make_hist([elt[1] for elt in mean_accy_copy if elt[0]=='healthy'],
             [elt[1] for elt in mean_accy_copy if elt[0]=='impaired'],
             binedges,'Mean y Acceleration','mean_accy_copy',path)
binedges = ct.get_bins(stddev_accx_copy,nbins=10)
ct.make_hist([elt[1] for elt in stddev_accx_copy if elt[0]=='healthy'],
             [elt[1] for elt in stddev_accx_copy if elt[0]=='impaired'],
             binedges,'Standard Deviation of x Acceleration','stddev_accx_copy',path)
binedges = ct.get_bins(stddev_accy_copy,nbins=10)
ct.make_hist([elt[1] for elt in stddev_accy_copy if elt[0]=='healthy'],
             [elt[1] for elt in stddev_accy_copy if elt[0]=='impaired'],
             binedges,'Standard Deviation of y Acceleration','stddev_accy_copy',path)
binedges = ct.get_bins(mean_pressure_copy,nbins=10)
ct.make_hist([elt[1] for elt in mean_pressure_copy if elt[0]=='healthy'],
             [elt[1] for elt in mean_pressure_copy if elt[0]=='impaired'],
             binedges,'Average Pressure','mean_pressure_copy',path)
binedges = ct.get_bins(stddev_pressure_copy,nbins=10)
ct.make_hist([elt[1] for elt in stddev_pressure_copy if elt[0]=='healthy'],
             [elt[1] for elt in stddev_pressure_copy if elt[0]=='impaired'],
             binedges,'Standard Deviation of Pressure','stddev_pressure_copy',path)
binedges = ct.get_bins(relsidelobe_x_copy,nbins=10)
ct.make_hist([elt[1] for elt in relsidelobe_x_copy if elt[0]=='healthy'],
             [elt[1] for elt in relsidelobe_x_copy if elt[0]=='impaired'],
             binedges,'Relative Side Lobe Amplitude [dB] for X[k]','relsidelobeamplitude_x_copy',path)
binedges = ct.get_bins(relsidelobe_y_copy,nbins=10)
ct.make_hist([elt[1] for elt in relsidelobe_y_copy if elt[0]=='healthy'],
             [elt[1] for elt in relsidelobe_y_copy if elt[0]=='impaired'],
             binedges,'Relative Side Lobe Amplitude [dB] for Y[k]','relsidelobeamplitude_y_copy',path)
# command clocks
binedges = ct.get_bins(mean_velx_command,nbins=10)
ct.make_hist([elt[1] for elt in mean_velx_command if elt[0]=='healthy'],
             [elt[1] for elt in mean_velx_command if elt[0]=='impaired'],
             binedges,'Mean x Velocity','mean_velx_command',path)
binedges = ct.get_bins(mean_vely_command,nbins=10)
ct.make_hist([elt[1] for elt in mean_vely_command if elt[0]=='healthy'],
             [elt[1] for elt in mean_vely_command if elt[0]=='impaired'],
             binedges,'Mean y Velocity','mean_vely_command',path)
binedges = ct.get_bins(stddev_velx_command,nbins=10)
ct.make_hist([elt[1] for elt in stddev_velx_command if elt[0]=='healthy'],
             [elt[1] for elt in stddev_velx_command if elt[0]=='impaired'],
             binedges,'Standard Deviation of x Velocity','stddev_velx_command',path)
binedges = ct.get_bins(stddev_vely_command,nbins=10)
ct.make_hist([elt[1] for elt in stddev_vely_command if elt[0]=='healthy'],
             [elt[1] for elt in stddev_vely_command if elt[0]=='impaired'],
             binedges,'Standard Deviation of y Velocity','stddev_vely_command',path)
binedges = ct.get_bins(mean_accx_command,nbins=10)
ct.make_hist([elt[1] for elt in mean_accx_command if elt[0]=='healthy'],
             [elt[1] for elt in mean_accx_command if elt[0]=='impaired'],
             binedges,'Mean x Acceleration','mean_accx_command',path)
binedges = ct.get_bins(mean_accy_command,nbins=10)
ct.make_hist([elt[1] for elt in mean_accy_command if elt[0]=='healthy'],
             [elt[1] for elt in mean_accy_command if elt[0]=='impaired'],
             binedges,'Mean y Acceleration','mean_accy_command',path)
binedges = ct.get_bins(stddev_accx_command,nbins=10)
ct.make_hist([elt[1] for elt in stddev_accx_command if elt[0]=='healthy'],
             [elt[1] for elt in stddev_accx_command if elt[0]=='impaired'],
             binedges,'Standard Deviation of x Acceleration','stddev_accx_command',path)
binedges = ct.get_bins(stddev_accy_command,nbins=10)
ct.make_hist([elt[1] for elt in stddev_accy_command if elt[0]=='healthy'],
             [elt[1] for elt in stddev_accy_command if elt[0]=='impaired'],
             binedges,'Standard Deviation of y Acceleration','stddev_accy_command',path)
binedges = ct.get_bins(mean_pressure_command,nbins=10)
ct.make_hist([elt[1] for elt in mean_pressure_command if elt[0]=='healthy'],
             [elt[1] for elt in mean_pressure_command if elt[0]=='impaired'],
             binedges,'Average Pressure','mean_pressure_command',path)
binedges = ct.get_bins(stddev_pressure_command,nbins=10)
ct.make_hist([elt[1] for elt in stddev_pressure_command if elt[0]=='healthy'],
             [elt[1] for elt in stddev_pressure_command if elt[0]=='impaired'],
             binedges,'Standard Deviation of Pressure','stddev_pressure_command',path)
binedges = ct.get_bins(relsidelobe_x_command,nbins=10)
ct.make_hist([elt[1] for elt in relsidelobe_x_command if elt[0]=='healthy'],
             [elt[1] for elt in relsidelobe_x_command if elt[0]=='impaired'],
             binedges,'Relative Side Lobe Amplitude [dB] for X[k]','relsidelobeamplitude_x_command',path)
binedges = ct.get_bins(relsidelobe_y_command,nbins=10)
ct.make_hist([elt[1] for elt in relsidelobe_y_command if elt[0]=='healthy'],
             [elt[1] for elt in relsidelobe_y_command if elt[0]=='impaired'],
             binedges,'Relative Side Lobe Amplitude [dB] for Y[k]','relsidelobeamplitude_y_command',path)
