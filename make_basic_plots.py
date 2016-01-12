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
        
    # get coordinates and timestamps
    x_copy,y_copy,p_copy,t_copy,x_command,y_command,p_command,t_command = ct.parse_file(path,fname)
    print '\n',len(x_copy),len(y_copy)
    print len(x_command),len(y_command),'\n'
    
    # interpolate to find missing points if necessary
    # assume for right now that we are never missing multiple consecutive points
    x_copy,y_copy = ct.interpolate(x_copy,y_copy)
    x_command,y_command = ct.interpolate(x_command,y_command)
    print '\n',len(x_copy),len(y_copy)
    print len(x_command),len(y_command),'\n'
    
    x_copy,y_copy,p_copy,t_copy = np.array(x_copy),np.array(y_copy),np.array(p_copy),np.array(t_copy)
    x_command,y_command,p_command,t_command = np.array(x_command),np.array(y_command),np.array(p_command),np.array(t_command)
    
    # plot x[t], y[t], v[t], and a[t]
    vx_copy,vy_copy,ax_copy,ay_copy = ct.deriv_doublederiv(x_copy,y_copy,t_copy)
    vx_command,vy_command,ax_command,ay_command = ct.deriv_doublederiv(x_command,y_command,t_command)
    ct.plot_xyt_other(x,y,t,vx_copy,vy_copy,t_copy,path,othername,fname)

    # plot spectrum of x[t] and y[t]
    
    # plot pressure
    
    print 'Done'
