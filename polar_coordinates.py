# -*- coding: utf-8 -*-
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

pi = math.pi

# copy or command clock?
clock_type = 'copy'

# path to polar coordinate figures
if not os.path.exists(path+'polar_coordinates'):
    os.makedirs(path+'polar_coordinates')
if not os.path.exists(path+'polar_coordinates/'+clock_type+'_healthy'):
    os.makedirs(path+'polar_coordinates/'+clock_type+'_healthy')
if not os.path.exists(path+'polar_coordinates/'+clock_type+'_impaired'):
    os.makedirs(path+'polar_coordinates/'+clock_type+'_impaired')

for fname in dirs:
    if 'CIN' not in fname and 'YDU' not in fname:
        continue
    print 'reading file ',fname,'...'
    
    ftype = ''
    if 'YDU' in fname:
        ftype = 'healthy'
    elif 'CIN' in fname:
        ftype = 'impaired'
    else:
        print 'not a valid file name'

    x_eqdist = np.loadtxt(path+'norm_velocity_data/'+fname[:len(fname)-4]+'_x_eqdist_'+clock_type+'.txt')
    y_eqdist = np.loadtxt(path+'norm_velocity_data/'+fname[:len(fname)-4]+'_y_eqdist_'+clock_type+'.txt')

    if len(x_eqdist)==0:
        continue
        
    # plot circle in polar coordinates
    # find COM
    x_com = np.mean(x_eqdist)
    y_com = np.mean(y_eqdist)

    # get r and theta
    r,theta = [],[]
    for w in range(len(x_eqdist)):
        dx,dy = x_eqdist[w]-x_com,y_eqdist[w]-y_com
        dist = sqrt(dx**2+dy**2)
        angle = math.atan2(dy,dx)
        r.append(dist)
        theta.append(angle)
    r,theta = np.array(r),np.array(theta)

    # figure settings
    x_axis_fontsize = 15
    y_axis_fontsize = 20

    # plot
    plt.close('all')
    fig_rtheta = plt.figure()
    rtheta = fig_rtheta.add_subplot(111)
    rtheta.plot(theta,r)

    # set axis limits
    rtheta.set_xlim(left=-3.2,right=3.2)

    # set axis labels
    rtheta.set_xlabel(r'$\theta$',fontsize=x_axis_fontsize)
    rtheta.set_ylabel(r'$r$',fontsize=y_axis_fontsize)

    # add drawing type (healthy or impaired) and file name
    fig_rtheta.text(0.99, 0.96,fname[:len(fname)-4],fontsize=10,color='red',va='baseline',ha='right',multialignment='left')

    if 'YDU' in fname:
        fig_rtheta.text(0.32, 0.955, 'HEALTHY ('+clock_type+')',fontsize=15,color='black',va='baseline',ha='right',multialignment='left')
    elif 'CIN' in fname:
        fig_rtheta.text(0.32, 0.955, 'IMPAIRED ('+clock_type+')',fontsize=15,color='black',va='baseline',ha='right',multialignment='left')
    else:
        print 'not a valid filename'

    # save figure
    fig_rtheta.savefig(path+'figs_raw/'+fname[:len(fname)-4]+'/xy_polar_'+clock_type+'_'+fname[:len(fname)-4]+'.png')
    fig_rtheta.savefig(path+'polar_coordinates/'+clock_type+'_'+ftype+'/xy_polar_'+clock_type+'_'+fname[:len(fname)-4]+'.png')
