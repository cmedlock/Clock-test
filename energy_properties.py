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

path = '/Users/cmedlock/Documents/DSP_UROP/DataCatherine/'
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

# estimation of eccentricity of ellipse
eccentricity = []

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

    f = open(path+fname)
    data = f.readlines()

    if not os.path.exists(path+'figs_raw/'+fname[:len(fname)-4]):
        os.makedirs(path+'figs_raw/'+fname[:len(fname)-4])

    x,y,t = [],[],[]

    # read in data
    record,found_clock = False,False
    for w in range(len(data)):
        line = data[w]
        # found copy clock?
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
    
    if len(x)==0:
        continue

    # change x so that if the whole clock is drawn,
    # it is oriented correctly
    x = max(x)+10-x
    
    # normalize the size of the word such that all coordinates
    # are linearly positioned between 0 and 127
    x = 127*(x-min(x))/(max(x)-min(x))
    y = 127*(y-min(y))/(max(y)-min(y))

    # compensate for non-constant velocity
    
    N_orig = len(x)
    N_new = 250

    # calculate average distance between points
    dists = []
    for w in range(1,len(x)):
        dx,dy = x[w]-x[w-1],y[w]-y[w-1]
        dist = math.sqrt(dx**2+dy**2)
        dists.append(dist)
    dist_avg = mean(dists)
    dist_total = sum(dists)
    #print 'average distance between points is ',dist_avg_copy
    #print 'total distance is ',sum(dists_copy)

    # now want to get N_orig evenly-spaced points along the curve

    # generate a much longer array with 199 linearly-interpolated 
    # points between the actual data points
    x_interp,y_interp = [],[]
    for w in range(len(x)-1):
        x_interp.append(x[w])
        y_interp.append(y[w])
        dx,dy = x[w+1]-x[w],y[w+1]-y[w]
        dist = math.sqrt(dx**2+dy**2)
        n_segments = ceil(dist/dist_avg)*200
        for r in range(1,int(n_segments)):
            x_new = x[w]+r*dx/n_segments
            y_new = y[w]+r*dy/n_segments
            x_interp.append(x_new)
            y_interp.append(y_new)
    x_interp.append(x[-1])
    y_interp.append(y[-1])

    # start from the first point and find the ones that are 
    # approximately a distance dist_avg from each other
    x_eqdist,y_eqdist = [x_interp[0]],[y_interp[0]]
    idx = 0
    for k in range(N_new):
        dist_sofar = 0
        for j in range(idx,len(x_interp)-1):
            dx,dy = x_interp[j+1]-x_interp[j],y_interp[j+1]-y_interp[j]
            dist_sofar += math.sqrt(dx**2+dy**2)
            #if abs(dist_sofar-dist_avg)<dist_avg/100.:
            if abs(dist_sofar-dist_total/250.)<dist_total/(250.*100.):
                idx = j+1
	        break
        x_eqdist.append(x_interp[idx])
        y_eqdist.append(y_interp[idx])

    # subtract mean values so there is no DC term adding an extra pole
    x_eqdist = [elt-mean(x_eqdist) for elt in x_eqdist]
    y_eqdist = [elt-mean(y_eqdist) for elt in y_eqdist]

    # dft's
    dft_size = len(x_eqdist)
    k = np.arange(dft_size)
    dftx,dfty = np.fft.fft(x_eqdist,n=dft_size),np.fft.fft(y_eqdist,n=dft_size)
    
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
    Ex,Ey = np.abs(dftx)**2,np.abs(dfty)**2
    Ex_total,Ey_total = sum(Ex),sum(Ey)
    Ex_peak,Ey_peak = 2*Ex[1]/Ex_total,2*Ey[1]/Ey_total
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
        #if angle<0:
            #angle = angle+2*pi
        r.append(dist)
        theta.append(angle)
    r,theta = np.array(r),np.array(theta)
    # estimate a,b of ellipse as avg values of 30 largest and 30 smallest r values
    r_sorted = r
    r_sorted.sort()
    a_estimate = mean(r_sorted[-10:])
    b_estimate = mean(r_sorted[:10])
    e_estimate = math.sqrt(1-(b_estimate/a_estimate)**2)
    eccentricity.append((ftype,e_estimate))

# compare energy properties for the drawings of healthy vs. impaired patients
#ct.make_hist([elt[1] for elt in Epeak_x if elt[0]=='healthy'],
#             [elt[1] for elt in Epeak_x if elt[0]=='impaired'],
#             10,'Fraction of Energy in Fundamental Freq. for x[n]','Epeak_x_'+clock_type,path)
#ct.make_hist([elt[1] for elt in Epeak_y if elt[0]=='healthy'],
#             [elt[1] for elt in Epeak_y if elt[0]=='impaired'],
#             10,'Fraction of Energy in Fundamental Freq. for y[n]','Epeak_y_'+clock_type,path)

ct.make_hist([elt[1] for elt in eccentricity if elt[0]=='healthy'],
             [elt[1] for elt in eccentricity if elt[0]=='impaired'],
             10,'Estimation of Ellipse Eccentricity','eccentricity_'+clock_type,path)

#ct.make_hist([elt[1] for elt in Estd_x if elt[0]=='healthy'],
#             [elt[1] for elt in Estd_x if elt[0]=='impaired'],
#             10,'Std. Deviation of Energy Distribution','Estd_x_'+clock_type,path)
#ct.make_hist([elt[1] for elt in Estd_y if elt[0]=='healthy'],
#             [elt[1] for elt in Estd_y if elt[0]=='impaired'],
#             10,'Std. Deviation of Energy Distribution','Estd_y_'+clock_type,path)
#
#ct.make_hist([elt[1] for elt in Ecentral_x if elt[0]=='healthy'],
#             [elt[1] for elt in Ecentral_x if elt[0]=='impaired'],
#             10,'Fraction of Energy w/in 1 Std.Dev. of w = 0','Ecentral_x_'+clock_type,path)
#ct.make_hist([elt[1] for elt in Ecentral_y if elt[0]=='healthy'],
#             [elt[1] for elt in Ecentral_y if elt[0]=='impaired'],
#             10,'Fraction of Energy w/in 1 Std.Dev. of w = 0','Ecentral_y_'+clock_type,path)
