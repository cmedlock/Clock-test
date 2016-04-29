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
clock_type = 'command'

# path to polar clock face figures
if not os.path.exists(path+'clock_face_drawings'):
    os.makedirs(path+'clock_face_drawings')
if not os.path.exists(path+'clock_face_drawings/'+clock_type+'_healthy'):
    os.makedirs(path+'clock_face_drawings/'+clock_type+'_healthy')
if not os.path.exists(path+'clock_face_drawings/'+clock_type+'_impaired'):
    os.makedirs(path+'clock_face_drawings/'+clock_type+'_impaired')

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
        
    # figure settings
    x_axis_fontsize = 15
    y_axis_fontsize = 20

    # plot
    fig_xy = plt.figure()
    xy = fig_xy.add_subplot(111)
    xy.plot(x_eqdist,y_eqdist,lw=2)

    # set axis labels
    xy.set_xlabel(r'$y$',fontsize=x_axis_fontsize)
    xy.set_ylabel(r'$x$',fontsize=y_axis_fontsize)

    # add drawing type (healthy or impaired) and file name
    fig_xy.text(0.99, 0.96,fname[:len(fname)-4],fontsize=10,color='red',va='baseline',ha='right',multialignment='left')

    if 'YDU' in fname:
        fig_xy.text(0.32, 0.955, 'HEALTHY ('+clock_type+')',fontsize=15,color='black',va='baseline',ha='right',multialignment='left')
    elif 'CIN' in fname:
        fig_xy.text(0.32, 0.955, 'IMPAIRED ('+clock_type+')',fontsize=15,color='black',va='baseline',ha='right',multialignment='left')
    else:
        print 'not a valid filename'

    # save figure
    fig_xy.savefig(path+'figs_raw/'+fname[:len(fname)-4]+'/xy_'+clock_type+'_'+fname[:len(fname)-4]+'.png')
    fig_xy.savefig(path+'clock_face_drawings/'+clock_type+'_'+ftype+'/xy_'+clock_type+'_'+fname[:len(fname)-4]+'.png')
