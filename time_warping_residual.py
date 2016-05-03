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

time_warping_resid_x,time_warping_resid_y = [],[]

# copy or command clock?
clock_type = 'command'

# path to polar coordinate figures
if not os.path.exists(path+'time_warping_residual'):
    os.makedirs(path+'time_warping_residual')
if not os.path.exists(path+'compare_healthy_impaired'):
    os.makedirs(path+'compare_healthy_impaired')

for fname in dirs:
    if 'CIN' not in fname and 'YDU' not in fname:
        continue
    #print 'reading file ',fname,'...'
    
    ftype = ''
    if 'YDU' in fname:
        ftype = 'healthy'
    elif 'CIN' in fname:
        ftype = 'impaired'
    else:
        print 'not a valid file name'

    f = open(path+fname)
    data = f.readlines()

    # get raw data (i.e. not time-warped)
    x,y,t = [],[],[]

    # read in data
    record,found_clock = False,False
    for w in range(len(data)):
        line = data[w]
        # found clock?
        if found_clock==False:
            if clock_type.upper() in line:
                found_clock = True
            continue
        # start recording?
        elif found_clock==True and 'CLOCKFACE' in line:
            record = True
            continue
        # stop recording?
        elif record==True:
            if 'symbol label' in line and len(x)>0:
                record = False
                break
            elif 'point' not in line:
                continue
        # done?
        elif found_clock==False and record==False and len(x)>0:
            break
        # other?
        else:
            continue
        line = line.split('"')
        
        # position and time
        xcoord = double(line[3])
        ycoord = double(line[1])
        timestamp = double(line[7])
        x.append(xcoord)
        y.append(ycoord)
        t.append(timestamp)
    
    f.close()

    x,y = np.array(x),np.array(y)
    
    # get time-warped data
    x_eqdist = np.loadtxt(path+'norm_velocity_data/'+fname[:len(fname)-4]+'_x_eqdist_'+clock_type+'.txt')
    y_eqdist = np.loadtxt(path+'norm_velocity_data/'+fname[:len(fname)-4]+'_y_eqdist_'+clock_type+'.txt')

    if len(x_eqdist)==0:
        continue
    
    # find energy in difference between raw data and time-warped data
    x_resid,y_resid = x-x_eqdist,y-y_eqdist
    time_warping_resid_x.append((ftype,sum(x_resid**2)*10**(-5)))
    time_warping_resid_y.append((ftype,sum(y_resid**2)*10**(-5)))

# compare energy in difference between raw data and time-warped data
# (twr = time-warping residual)
twr_x_all = [elt[1] for elt in time_warping_resid_x]
mean_x,std_x = mean(twr_x_all),std(twr_x_all)
twr_y_all = [elt[1] for elt in time_warping_resid_y]
mean_y,std_y = mean(twr_y_all),std(twr_y_all)
ct.make_hist([elt[1] for elt in time_warping_resid_x if elt[0]=='healthy'],
             [elt[1] for elt in time_warping_resid_x if elt[0]=='impaired'],
             10,mean_x-std_x,mean_x+std_x,r'Energy in Time-Warped Residual $\times 10^{-5}$ for x[n]','time_warping_resid_energy_x_'+clock_type,path)
ct.make_hist([elt[1] for elt in time_warping_resid_y if elt[0]=='healthy'],
             [elt[1] for elt in time_warping_resid_y if elt[0]=='impaired'],
             10,mean_y-std_y,mean_y+2*std_y,r'Energy in Time-Warped Residual $\times 10^{-5}$ for y[n]','time_warping_resid_energy_y_'+clock_type,path)
# in case the histograms don't come out right
np.savetxt(path+'time_warping_residual/time_warping_resid_energy_x_healthy_'+clock_type+'.txt',[elt[1] for elt in time_warping_resid_x if elt[0]=='healthy'])
np.savetxt(path+'time_warping_residual/time_warping_resid_energy_x_impaired_'+clock_type+'.txt',[elt[1] for elt in time_warping_resid_x if elt[0]=='impaired'])
np.savetxt(path+'time_warping_residual/time_warping_resid_energy_y_healthy_'+clock_type+'.txt',[elt[1] for elt in time_warping_resid_y if elt[0]=='healthy'])
np.savetxt(path+'time_warping_residual/time_warping_resid_energy_y_impaired_'+clock_type+'.txt',[elt[1] for elt in time_warping_resid_y if elt[0]=='impaired'])
        