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

path = '/Users/cmedlock/Documents/DSP_UROP/DataCatherine/'
dirs = os.listdir(path)

pi = math.pi

if not os.path.exists(path+'norm_velocity_data'):
    os.makedirs(path+'norm_velocity_data')

# copy or command clock?
clock_type = 'COPY'

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

    x,y,t = [],[],[]

    # read in data
    record,found_clock = False,False
    for w in range(len(data)):
        line = data[w]
        # found clock?
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
        np.savetxt(path+'norm_velocity_data/'+fname[:len(fname)-4]+'_x_eqdist_'+clock_type+'.txt',x)
        np.savetxt(path+'norm_velocity_data/'+fname[:len(fname)-4]+'_y_eqdist_'+clock_type+'.txt',y)
        continue

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

    # if the points are already evenly-spaced, don't interpolate
    if np.var(np.array(dists))<10**-12:
        x_eqdist,y_eqdist = x,y
    else:
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
        for k in range(N_new-1):
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
    np.savetxt(path+'norm_velocity_data/'+fname[:len(fname)-4]+'_x_eqdist_'+clock_type+'.txt',x_eqdist)
    np.savetxt(path+'norm_velocity_data/'+fname[:len(fname)-4]+'_y_eqdist_'+clock_type+'.txt',y_eqdist)
