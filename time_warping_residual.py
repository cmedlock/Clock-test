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

# (twr = time-warping residual)
time_warping_resid_x,time_warping_resid_y = [],[]
twr_q1_x,twr_q1_y = [],[]
twr_q2_x,twr_q2_y = [],[]
twr_q3_x,twr_q3_y = [],[]
twr_q4_x,twr_q4_y = [],[]

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
    
    # separate the energy into different quadrants
    Eresid_q1_x,Eresid_q2_x,Eresid_q3_x,Eresid_q4_x = 0,0,0,0
    Eresid_q1_y,Eresid_q2_y,Eresid_q3_y,Eresid_q4_y = 0,0,0,0
    for w in range(len(x_resid)):
        if x_eqdist[w]>0 and y_eqdist[w]>0:
            Eresid_q1_x += x_resid[w]**2*10**(-5)
            Eresid_q1_y += y_resid[w]**2*10**(-5)
        elif x_eqdist[w]<0 and y_eqdist[w]>0:
            Eresid_q2_x += x_resid[w]**2*10**(-5)
            Eresid_q2_y += y_resid[w]**2*10**(-5)
        elif x_eqdist[w]<0 and y_eqdist[w]<0:
            Eresid_q3_x += x_resid[w]**2*10**(-5)
            Eresid_q3_y += y_resid[w]**2*10**(-5)
        elif x_eqdist[w]>0 and y_eqdist[w]<0:
            Eresid_q4_x += x_resid[w]**2*10**(-5)
            Eresid_q4_y += y_resid[w]**2*10**(-5)
    twr_q1_x.append((ftype,Eresid_q1_x))
    twr_q1_y.append((ftype,Eresid_q1_y))
    twr_q2_x.append((ftype,Eresid_q2_x))
    twr_q2_y.append((ftype,Eresid_q2_y))
    twr_q3_x.append((ftype,Eresid_q3_x))
    twr_q3_y.append((ftype,Eresid_q3_y))
    twr_q4_x.append((ftype,Eresid_q4_x))
    twr_q4_y.append((ftype,Eresid_q4_y))
            

# compare energy in difference between raw data and time-warped data
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
# quadrant 1
twr_q1_x_all = [elt[1] for elt in twr_q1_x]
mean_q1_x,std_q1_x = mean(twr_q1_x_all),std(twr_q1_x_all)
twr_q1_y_all = [elt[1] for elt in twr_q1_y]
mean_q1_y,std_q1_y = mean(twr_q1_y_all),std(twr_q1_y_all)
ct.make_hist([elt[1] for elt in twr_q1_x if elt[0]=='healthy'],
             [elt[1] for elt in twr_q1_x if elt[0]=='impaired'],
             10,mean_q1_x-std_q1_x,mean_q1_x+std_q1_x,r'Quadrant 1: Energy in Time-Warped Residual $\times 10^{-5}$ for x[n]','twr_q1_energy_x_'+clock_type,path)
ct.make_hist([elt[1] for elt in twr_q1_y if elt[0]=='healthy'],
             [elt[1] for elt in twr_q1_y if elt[0]=='impaired'],
             10,mean_q1_y-std_q1_y,mean_q1_y+std_q1_y,r'Quadrant 1: Energy in Time-Warped Residual $\times 10^{-5}$ for y[n]','twr_q1_energy_y_'+clock_type,path)
# in case the histograms don't come out right
np.savetxt(path+'time_warping_residual/twr_q1_energy_x_healthy_'+clock_type+'.txt',[elt[1] for elt in twr_q1_x if elt[0]=='healthy'])
np.savetxt(path+'time_warping_residual/twr_q1_energy_x_impaired_'+clock_type+'.txt',[elt[1] for elt in twr_q1_x if elt[0]=='impaired'])
np.savetxt(path+'time_warping_residual/twr_q1_energy_y_healthy_'+clock_type+'.txt',[elt[1] for elt in twr_q1_y if elt[0]=='healthy'])
np.savetxt(path+'time_warping_residual/twr_q1_energy_y_impaired_'+clock_type+'.txt',[elt[1] for elt in twr_q1_y if elt[0]=='impaired'])
# quadrant 2
twr_q2_x_all = [elt[1] for elt in twr_q2_x]
mean_q2_x,std_q2_x = mean(twr_q2_x_all),std(twr_q2_x_all)
twr_q2_y_all = [elt[1] for elt in twr_q2_y]
mean_q2_y,std_q2_y = mean(twr_q2_y_all),std(twr_q2_y_all)
ct.make_hist([elt[1] for elt in twr_q2_x if elt[0]=='healthy'],
             [elt[1] for elt in twr_q2_x if elt[0]=='impaired'],
             10,mean_q2_x-std_q2_x,mean_q2_x+std_q2_x,r'Quadrant 1: Energy in Time-Warped Residual $\times 10^{-5}$ for x[n]','twr_q2_energy_x_'+clock_type,path)
ct.make_hist([elt[1] for elt in twr_q2_y if elt[0]=='healthy'],
             [elt[1] for elt in twr_q2_y if elt[0]=='impaired'],
             10,mean_q2_y-std_q2_y,mean_q2_y+std_q2_y,r'Quadrant 1: Energy in Time-Warped Residual $\times 10^{-5}$ for y[n]','twr_q2_energy_y_'+clock_type,path)
# in case the histograms don't come out right
np.savetxt(path+'time_warping_residual/twr_q2_energy_x_healthy_'+clock_type+'.txt',[elt[1] for elt in twr_q2_x if elt[0]=='healthy'])
np.savetxt(path+'time_warping_residual/twr_q2_energy_x_impaired_'+clock_type+'.txt',[elt[1] for elt in twr_q2_x if elt[0]=='impaired'])
np.savetxt(path+'time_warping_residual/twr_q2_energy_y_healthy_'+clock_type+'.txt',[elt[1] for elt in twr_q2_y if elt[0]=='healthy'])
np.savetxt(path+'time_warping_residual/twr_q2_energy_y_impaired_'+clock_type+'.txt',[elt[1] for elt in twr_q2_y if elt[0]=='impaired'])
# quadrant 3
twr_q3_x_all = [elt[1] for elt in twr_q3_x]
mean_q3_x,std_q3_x = mean(twr_q3_x_all),std(twr_q3_x_all)
twr_q3_y_all = [elt[1] for elt in twr_q3_y]
mean_q3_y,std_q3_y = mean(twr_q3_y_all),std(twr_q3_y_all)
ct.make_hist([elt[1] for elt in twr_q3_x if elt[0]=='healthy'],
             [elt[1] for elt in twr_q3_x if elt[0]=='impaired'],
             10,mean_q3_x-std_q3_x,mean_q3_x+std_q3_x,r'Quadrant 1: Energy in Time-Warped Residual $\times 10^{-5}$ for x[n]','twr_q3_energy_x_'+clock_type,path)
ct.make_hist([elt[1] for elt in twr_q3_y if elt[0]=='healthy'],
             [elt[1] for elt in twr_q3_y if elt[0]=='impaired'],
             10,mean_q3_y-std_q3_y,mean_q3_y+std_q3_y,r'Quadrant 1: Energy in Time-Warped Residual $\times 10^{-5}$ for y[n]','twr_q3_energy_y_'+clock_type,path)
# in case the histograms don't come out right
np.savetxt(path+'time_warping_residual/twr_q3_energy_x_healthy_'+clock_type+'.txt',[elt[1] for elt in twr_q3_x if elt[0]=='healthy'])
np.savetxt(path+'time_warping_residual/twr_q3_energy_x_impaired_'+clock_type+'.txt',[elt[1] for elt in twr_q3_x if elt[0]=='impaired'])
np.savetxt(path+'time_warping_residual/twr_q3_energy_y_healthy_'+clock_type+'.txt',[elt[1] for elt in twr_q3_y if elt[0]=='healthy'])
np.savetxt(path+'time_warping_residual/twr_q3_energy_y_impaired_'+clock_type+'.txt',[elt[1] for elt in twr_q3_y if elt[0]=='impaired'])
# quadrant 4
twr_q4_x_all = [elt[1] for elt in twr_q4_x]
mean_q4_x,std_q4_x = mean(twr_q4_x_all),std(twr_q4_x_all)
twr_q4_y_all = [elt[1] for elt in twr_q4_y]
mean_q4_y,std_q4_y = mean(twr_q4_y_all),std(twr_q4_y_all)
ct.make_hist([elt[1] for elt in twr_q4_x if elt[0]=='healthy'],
             [elt[1] for elt in twr_q4_x if elt[0]=='impaired'],
             10,mean_q4_x-std_q4_x,mean_q4_x+std_q4_x,r'Quadrant 1: Energy in Time-Warped Residual $\times 10^{-5}$ for x[n]','twr_q4_energy_x_'+clock_type,path)
ct.make_hist([elt[1] for elt in twr_q4_y if elt[0]=='healthy'],
             [elt[1] for elt in twr_q4_y if elt[0]=='impaired'],
             10,mean_q4_y-std_q4_y,mean_q4_y+std_q4_y,r'Quadrant 1: Energy in Time-Warped Residual $\times 10^{-5}$ for y[n]','twr_q4_energy_y_'+clock_type,path)
# in case the histograms don't come out right
np.savetxt(path+'time_warping_residual/twr_q4_energy_x_healthy_'+clock_type+'.txt',[elt[1] for elt in twr_q4_x if elt[0]=='healthy'])
np.savetxt(path+'time_warping_residual/twr_q4_energy_x_impaired_'+clock_type+'.txt',[elt[1] for elt in twr_q4_x if elt[0]=='impaired'])
np.savetxt(path+'time_warping_residual/twr_q4_energy_y_healthy_'+clock_type+'.txt',[elt[1] for elt in twr_q4_y if elt[0]=='healthy'])
np.savetxt(path+'time_warping_residual/twr_q4_energy_y_impaired_'+clock_type+'.txt',[elt[1] for elt in twr_q4_y if elt[0]=='impaired'])      
