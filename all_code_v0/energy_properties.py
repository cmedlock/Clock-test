# calculate DFS coefficients of each x vs. t and y vs. t signal
# calculate percent energy in peak, std. deviation of energy distribution,
# and percent energy within 1 std. deviation of w = 0

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

# copy or command clock?
clock_type = 'COPY'

# save interesting quantities
# fraction of energy contained in largest DFS coefficient
Epeak_x,Epeak_y = [],[]

# fraction of energy contained in DFS coefficients that
# within 1 std. deviation of w = 0
Estd_x,Estd_y = [],[]
Ecentral_x,Ecentral_y = [],[]

for fname in dirs:
    if 'Scored' not in fname:
        continue
    print 'reading file ',fname,'...'
    
    ftype = ''
    if 'YDU' in fname:
        ftype = 'healthy'
    elif 'CIN' in fname:
        ftype = 'impaired'
    else:
        print 'not a valid file name'

    f = open(path+fname)
    data = f.readlines()

    if not os.path.exists(path+'figs_raw/'+fname[:len(fname)-4]):
        os.makedirs(path+'figs_raw/'+fname[:len(fname)-4])

    x,y,t = [],[],[]

    # read in data
    record,found_clock = False,False
    for w in range(len(data)):
        line = data[w]
        # found lock?
        if found_clock==False:
            if clock_type in line:
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

    x,y,t = np.array(x),np.array(y),np.array(t)

    # subtract mean values (zm = zero mean)
    x_zm,y_zm = x-mean(x),y-mean(y)
    
    # dft's
    dft_size = len(x)
    k = np.arange(dft_size)
    dftx,dfty = np.fft.fft(x_zm,n=dft_size),np.fft.fft(y_zm,n=dft_size)
    
    # k_near_pi is the smallest k value for which w_k = 2*pi*k/N is
    # greater than or equal to pi
    k_near_pi = 0
    if dft_size%2==0:
        k_near_pi = dft_size/2+1
    else:
        k_near_pi = math.ceil(dft_size/2)
    
    # center the coefficients around w = 0
    k_centered = np.linspace(-dft_size/2,dft_size/2,dft_size)
    dftx_centered = np.concatenate((dftx[k_near_pi:],dftx[:k_near_pi]))
    dfty_centered = np.concatenate((dfty[k_near_pi:],dfty[:k_near_pi]))
    
    # percent energy in peak
    Ex,Ey = np.abs(dftx_centered)**2,np.abs(dfty_centered)**2
    Ex_total,Ey_total = sum(Ex),sum(Ey)
    Ex_peak,Ey_peak = 2*max(Ex)/Ex_total,2*max(Ey)/Ey_total
    Epeak_x.append((ftype,Ex_peak))
    Epeak_y.append((ftype,Ey_peak))

    # percent energy within 1 std. deviation of the center of the energy distribution
    Ex_var,Ey_var = 0,0
    for w in range(dft_size):
        Ex_var += k_centered[w]**2*Ex[w]/Ex_total
        Ey_var += k_centered[w]**2*Ey[w]/Ey_total
    Ex_std,Ey_std = math.sqrt(Ex_var),math.sqrt(Ey_var)
    Estd_x.append((ftype,Ex_std))
    Estd_y.append((ftype,Ey_std))
    
    Ex_central,Ey_central = 0,0
    for d in range(dft_size):
        if abs((dft_size-1)/2-d)<=Ex_std:
            Ex_central += Ex[d]/Ex_total
        if abs((dft_size-1)/2-d)<=Ey_std:
            Ey_central += Ey[d]/Ey_total
    Ecentral_x.append((ftype,Ex_central))
    Ecentral_y.append((ftype,Ey_central))

# compare energy properties for the drawings of healthy vs. impaired patients
ct.make_hist([elt[1] for elt in Epeak_x if elt[0]=='healthy'],
             [elt[1] for elt in Epeak_x if elt[0]=='impaired'],
             10,'Fraction of Energy in Largest DFS Coefficient','Epeak_x_'+clock_type,path)
ct.make_hist([elt[1] for elt in Epeak_y if elt[0]=='healthy'],
             [elt[1] for elt in Epeak_y if elt[0]=='impaired'],
             10,'Fraction of Energy in Largest DFS Coefficient','Epeak_y_'+clock_type,path)

ct.make_hist([elt[1] for elt in Estd_x if elt[0]=='healthy'],
             [elt[1] for elt in Estd_x if elt[0]=='impaired'],
             10,'Std. Deviation of Energy Distribution','Estd_x_'+clock_type,path)
ct.make_hist([elt[1] for elt in Estd_y if elt[0]=='healthy'],
             [elt[1] for elt in Estd_y if elt[0]=='impaired'],
             10,'Std. Deviation of Energy Distribution','Estd_y_'+clock_type,path)

ct.make_hist([elt[1] for elt in Ecentral_x if elt[0]=='healthy'],
             [elt[1] for elt in Ecentral_x if elt[0]=='impaired'],
             10,'Fraction of Energy w/in 1 Std.Dev. of w = 0','Ecentral_x_'+clock_type,path)
ct.make_hist([elt[1] for elt in Ecentral_y if elt[0]=='healthy'],
             [elt[1] for elt in Ecentral_y if elt[0]=='impaired'],
             10,'Fraction of Energy w/in 1 Std.Dev. of w = 0','Ecentral_y_'+clock_type,path)
