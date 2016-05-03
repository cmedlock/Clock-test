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
clock_type = 'command'

plt.close('all')
fig = plt.figure()
hist = fig.add_subplot(111)
nbins = 15

# percent energy in fundamental frequency for x[n]
Epeak_x_healthy = np.loadtxt(path+'DFS_coefficients/Epeak_x_healthy_'+clock_type+'.txt')
Epeak_x_impaired = np.loadtxt(path+'DFS_coefficients/Epeak_x_impaired_'+clock_type+'.txt')
Epeak_x_all = np.concatenate((Epeak_x_healthy,Epeak_x_impaired))
range_min,range_max = 0.97,1.001
ct.make_hist(Epeak_x_healthy,Epeak_x_impaired,15,range_min,range_max,'Relative Energy in Fundamental Freq. for x[n]','Epeak_x_'+clock_type,path)

# percent energy in fundamental frequency for y[n]
Epeak_y_healthy = np.loadtxt(path+'DFS_coefficients/Epeak_y_healthy_'+clock_type+'.txt')
Epeak_y_impaired = np.loadtxt(path+'DFS_coefficients/Epeak_y_impaired_'+clock_type+'.txt')
Epeak_y_all = np.concatenate((Epeak_y_healthy,Epeak_y_impaired))
range_min,range_max = 0.96,1.001
ct.make_hist(Epeak_y_healthy,Epeak_y_impaired,15,range_min,range_max,'Relative Energy in Fundamental Freq. for y[n]','Epeak_y_'+clock_type,path)

# percent energy in linear prediction error for x[n]
Eg_x_healthy = np.loadtxt(path+'lpc/Eg_x_healthy_'+clock_type+'.txt')
Eg_x_impaired = np.loadtxt(path+'lpc/Eg_x_impaired_'+clock_type+'.txt')
Eg_x_all = np.concatenate((Eg_x_healthy,Eg_x_impaired))
#range_min,range_max = 0,0.15 # copy clock
range_min,range_max = 0,0.2 # command clock
ct.make_hist(Eg_x_healthy,Eg_x_impaired,15,range_min,range_max,'Relative Energy in Linear Prediction Error for x[n]','Eg_x_'+clock_type,path)

# percent energy in linear prediction error for y[n]
Eg_y_healthy = np.loadtxt(path+'lpc/Eg_y_healthy_'+clock_type+'.txt')
Eg_y_impaired = np.loadtxt(path+'lpc/Eg_y_impaired_'+clock_type+'.txt')
Eg_y_all = np.concatenate((Eg_y_healthy,Eg_y_impaired))
range_min,range_max = 0,0.40
ct.make_hist(Eg_y_healthy,Eg_y_impaired,15,range_min,range_max,'Relative Energy in Linear Prediction Error for y[n]','Eg_y_'+clock_type,path)

# percent energy in complex linear prediction error for x[n]
Eg_x_healthy = np.loadtxt(path+'lpc_complex/Eg_x_healthy_'+clock_type+'.txt')
Eg_x_impaired = np.loadtxt(path+'lpc_complex/Eg_x_impaired_'+clock_type+'.txt')
Eg_x_all = np.concatenate((Eg_x_healthy,Eg_x_impaired))
#range_min,range_max = 0,0.15 # copy clock
range_min,range_max = 0,0.15 # command clock
ct.make_hist(Eg_x_healthy,Eg_x_impaired,15,range_min,range_max,'Relative Energy in Linear Prediction Error for x[n]','complex_Eg_x_'+clock_type,path)

# percent energy in complex linear prediction error for y[n]
Eg_y_healthy = np.loadtxt(path+'lpc_complex/Eg_y_healthy_'+clock_type+'.txt')
Eg_y_impaired = np.loadtxt(path+'lpc_complex/Eg_y_impaired_'+clock_type+'.txt')
Eg_y_all = np.concatenate((Eg_y_healthy,Eg_y_impaired))
range_min,range_max = 0,0.40
ct.make_hist(Eg_y_healthy,Eg_y_impaired,15,range_min,range_max,'Relative Energy in Linear Prediction Error for y[n]','complex_Eg_y_'+clock_type,path)

# short time frequency for x[n]
freq_var_x_healthy = np.loadtxt(path+'short_time_lp/freq_var_x_healthy_'+clock_type+'.txt')
freq_var_x_impaired = np.loadtxt(path+'short_time_lp/freq_var_x_impaired_'+clock_type+'.txt')
freq_var_x_all = np.concatenate((freq_var_x_healthy,freq_var_x_impaired))
#range_min,range_max = 0,0.01 # copy clock
range_min,range_max = 0,0.015 # command clock
ct.make_hist(freq_var_x_healthy,freq_var_x_impaired,15,range_min,range_max,'Variance in Short Time Frequency for x[n]','freq_var_x_'+clock_type,path)

# short time frequency for y[n]
freq_var_y_healthy = np.loadtxt(path+'short_time_lp/freq_var_y_healthy_'+clock_type+'.txt')
freq_var_y_impaired = np.loadtxt(path+'short_time_lp/freq_var_y_impaired_'+clock_type+'.txt')
freq_var_y_all = np.concatenate((freq_var_y_healthy,freq_var_y_impaired))
#range_min,range_max = 0,3500 # copy clock
range_min,range_max = 0,4000 # command clock
ct.make_hist(freq_var_y_healthy,freq_var_y_impaired,15,range_min,range_max,'Variance in Short Time Frequency for y[n]','freq_var_y_'+clock_type,path)
