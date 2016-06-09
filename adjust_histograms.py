# -*- coding: utf-8 -*-
import math
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import numpy.linalg
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
print 'clock_type is ',clock_type

plt.close('all')
fig = plt.figure()
hist = fig.add_subplot(111)
nbins = 15
"""
# a_1
print 'a_1'
a1_x_healthy = np.loadtxt(path+'lpc/a1_x_healthy_'+clock_type+'.txt')
a1_x_impaired = np.loadtxt(path+'lpc/a1_x_impaired_'+clock_type+'.txt')
a1_x_all = np.concatenate((a1_x_healthy,a1_x_impaired))
range_min_x,range_max_x = 1.6,2.05
a1_y_healthy = np.loadtxt(path+'lpc/a1_y_healthy_'+clock_type+'.txt')
a1_y_impaired = np.loadtxt(path+'lpc/a1_y_impaired_'+clock_type+'.txt')
a1_y_all = np.concatenate((a1_y_healthy,a1_y_impaired))
range_min_y,range_max_y = 1.3,2.1
ct.make_hist_xy(a1_x_healthy,a1_x_impaired,a1_y_healthy,a1_y_impaired,
                15,range_min_x,range_max_x,range_min_y,range_max_y,'a_1','a1_'+clock_type,path)

# a_2
print 'a_2'
a2_x_healthy = np.loadtxt(path+'lpc/a2_x_healthy_'+clock_type+'.txt')
a2_x_impaired = np.loadtxt(path+'lpc/a2_x_impaired_'+clock_type+'.txt')
a2_x_all = np.concatenate((a2_x_healthy,a2_x_impaired))
range_min_x,range_max_x = -1,-0.55
a2_y_healthy = np.loadtxt(path+'lpc/a2_y_healthy_'+clock_type+'.txt')
a2_y_impaired = np.loadtxt(path+'lpc/a2_y_impaired_'+clock_type+'.txt')
a2_y_all = np.concatenate((a2_y_healthy,a2_y_impaired))
range_min_y,range_max_y = -1,0
ct.make_hist_xy(a2_x_healthy,a2_x_impaired,a2_y_healthy,a2_y_impaired,
                15,range_min_x,range_max_x,range_min_y,range_max_y,'a_2','a2_'+clock_type,path)

# percent energy in linear prediction with timing error
print 'percent energy in linear prediction with timing error'
Eg_x_healthy = np.loadtxt(path+'lpc_complex_with_timing/Eg_x_healthy_'+clock_type+'.txt')
Eg_x_impaired = np.loadtxt(path+'lpc_complex_with_timing/Eg_x_impaired_'+clock_type+'.txt')
Eg_x_all = np.concatenate((Eg_x_healthy,Eg_x_impaired))
#range_min_x,range_max_x = 0,0.3 # command clock
range_min_x,range_max_x = 0,0.06 # copy clock
Eg_y_healthy = np.loadtxt(path+'lpc_complex_with_timing/Eg_y_healthy_'+clock_type+'.txt')
Eg_y_impaired = np.loadtxt(path+'lpc_complex_with_timing/Eg_y_impaired_'+clock_type+'.txt')
Eg_y_all = np.concatenate((Eg_y_healthy,Eg_y_impaired))
#range_min_y,range_max_y = 0,0.03 # command clock
range_min_y,range_max_y = 0,0.4 # copy clock
ct.make_hist_xy(Eg_x_healthy,Eg_x_impaired,Eg_y_healthy,Eg_y_impaired,
                15,range_min_x,range_max_x,range_min_y,range_max_y,'E_g','complex_Eg_with_timing_'+clock_type,path)

# percent energy in fundamental frequency
print 'percent energy in fundamental frequency:'
fig.clear()
Epeak_x_healthy = np.loadtxt(path+'DFS_coefficients/Epeak_x_healthy_'+clock_type+'.txt')
Epeak_x_impaired = np.loadtxt(path+'DFS_coefficients/Epeak_x_impaired_'+clock_type+'.txt')
Epeak_x_all = np.concatenate((Epeak_x_healthy,Epeak_x_impaired))
range_min_x,range_max_x = 0.96,1.001
Epeak_y_healthy = np.loadtxt(path+'DFS_coefficients/Epeak_y_healthy_'+clock_type+'.txt')
Epeak_y_impaired = np.loadtxt(path+'DFS_coefficients/Epeak_y_impaired_'+clock_type+'.txt')
Epeak_y_all = np.concatenate((Epeak_y_healthy,Epeak_y_impaired))
range_min_y,range_max_y = 0.96,1.001
ct.make_hist_xy(Epeak_x_healthy,Epeak_x_impaired,Epeak_y_healthy,Epeak_y_impaired,
                15,range_min_x,range_max_x,range_min_y,range_max_y,'E_fund','Epeak_'+clock_type,path)

# percent energy in linear prediction error
print 'percent energy in linear prediction error'
fig.clear()
Eg_x_healthy = np.loadtxt(path+'lpc/Eg_x_healthy_'+clock_type+'.txt')
Eg_x_impaired = np.loadtxt(path+'lpc/Eg_x_impaired_'+clock_type+'.txt')
Eg_x_all = np.concatenate((Eg_x_healthy,Eg_x_impaired))
range_min_x,range_max_x = 0,0.17
Eg_y_healthy = np.loadtxt(path+'lpc/Eg_y_healthy_'+clock_type+'.txt')
Eg_y_impaired = np.loadtxt(path+'lpc/Eg_y_impaired_'+clock_type+'.txt')
Eg_y_all = np.concatenate((Eg_y_healthy,Eg_y_impaired))
range_min_y,range_max_y = 0,0.40
ct.make_hist_xy(Eg_x_healthy,Eg_x_impaired,Eg_y_healthy,Eg_y_impaired,
            15,range_min_x,range_max_x,range_min_y,range_max_y,'E_g','Eg_'+clock_type,path)

# percent energy in complex linear prediction error
print 'percent energy in complex linear prediction error'
fig.clear()
Eg_x_healthy = np.loadtxt(path+'lpc_complex/Eg_x_healthy_'+clock_type+'.txt')
Eg_x_impaired = np.loadtxt(path+'lpc_complex/Eg_x_impaired_'+clock_type+'.txt')
Eg_x_all = np.concatenate((Eg_x_healthy,Eg_x_impaired))
range_min_x,range_max_x = 0,0.15 # copy clock
#range_min_x,range_max_x = 0,0.15 # command clock
Eg_y_healthy = np.loadtxt(path+'lpc_complex/Eg_y_healthy_'+clock_type+'.txt')
Eg_y_impaired = np.loadtxt(path+'lpc_complex/Eg_y_impaired_'+clock_type+'.txt')
Eg_y_all = np.concatenate((Eg_y_healthy,Eg_y_impaired))
range_min_y,range_max_y = 0,0.4
ct.make_hist_xy(Eg_x_healthy,Eg_x_impaired,Eg_y_healthy,Eg_y_impaired,
                15,range_min_x,range_max_x,range_min_y,range_max_y,'E_g','complex_Eg_'+clock_type,path)
"""
# short time frequency
print 'short time frequency'
fig.clear()
freq_var_x_healthy = np.loadtxt(path+'short_time_lp/freq_var_x_healthy_'+clock_type+'.txt')
freq_var_x_impaired = np.loadtxt(path+'short_time_lp/freq_var_x_impaired_'+clock_type+'.txt')
freq_var_x_all = np.concatenate((freq_var_x_healthy,freq_var_x_impaired))
range_min_x,range_max_x = 0,0.001 # copy clock
#range_min_x,range_max_x = 0,0.0005 # command clock
freq_var_y_healthy = np.loadtxt(path+'short_time_lp/freq_var_y_healthy_'+clock_type+'.txt')
freq_var_y_impaired = np.loadtxt(path+'short_time_lp/freq_var_y_impaired_'+clock_type+'.txt')
freq_var_y_all = np.concatenate((freq_var_y_healthy,freq_var_y_impaired))
range_min_y,range_max_y = 0,0.09
ct.make_hist_xy(freq_var_x_healthy,freq_var_x_impaired,freq_var_y_healthy,freq_var_y_impaired,
                15,range_min_x,range_max_x,range_min_y,range_max_y,'var(f(n))','freq_var_'+clock_type,path)