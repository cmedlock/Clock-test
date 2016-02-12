# -*- coding: utf-8 -*-
# make a video of the clock being drawn, from an angle different than the
# most straightforward one (looking straight down the time axis at the xy
# plane
import math
import matplotlib
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
import numpy as np
import os
from pylab import *

# path to data
path = '/Users/cmedlock/Documents/DSP_UROP/all_data/'
dirs = os.listdir(path)

# path to output
if not os.path.exists(path+'figs_raw'):
    os.makedirs(path+'figs_raw')

# make figs
for fname in dirs[:2]:
    if 'Scored' not in fname:
        continue
    print 'reading file ',fname,'...'
    f = open(path+fname)
    data = f.readlines()

    if not os.path.exists(path+'figs_raw/'+fname[:len(fname)-4]):
        os.makedirs(path+'figs_raw/'+fname[:len(fname)-4])

    # copy or command clock?
    clock_type = 'COMMAND'
    x,y,t = [],[],[]
    x_temp,y_temp,t_temp = [],[],[]

    # read in data
    record,found_clock,is_noise = False,False,False
    for w in range(len(data)):
        line = data[w]
        # found pen stroke?
        if found_clock==False and clock_type in line:
            found_clock = True
            continue
        # start recording?
        elif found_clock==True and 'CLOCKFACE' in line:
            record = True
            continue
        # stop recording?
        elif record==True:
            if 'RawData' in line and len(x)>0:
                break
            if 'symbol label' in line and len(x_temp)==0 and 'NOISE' in line:
                is_noise = True
            if 'symbol label' in line and len(x_temp)>0: # length requirement since sometimes there are empty symbols
                #print 'is_noise is now ',is_noise
                #print 'len(t_temp) is now ',len(t_temp)
                #print 'next symbol is ',line
                # store previous symbol and reset lists for the next one
                if is_noise==False:
                    #print 'storing previous symbol from ',t_temp[0],' to ',t_temp[-1],' (',(t_temp[0]-1307115563146.0)/10.,' to ',(t_temp[-1]-1307115563146.0)/10.,')'
                    x.append(x_temp)
                    y.append(y_temp)
                    t.append(t_temp)
                elif is_noise==True:
                    is_noise = False
                #print 'len(t) is now ',len(t)
                x_temp,y_temp,t_temp = [],[],[]
                if 'NOISE' in line:
                    is_noise = True
                continue
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
        x_temp.append(xcoord)
        y_temp.append(ycoord)
        t_temp.append(timestamp)
    
    f.close()

    # force all coordinates to be non-negative
    xmin = min([min(elt) for elt in x])
    ymin = min([min(elt) for elt in y])
    tmin = min([min(elt) for elt in t])
    x = [list(np.array(elt)-xmin) for elt in x]
    y = [list(np.array(elt)-ymin) for elt in y]
    t = [list((np.array(elt)-tmin)/10.) for elt in t] # the time array should have reasonable values
    
    # make each of x, y, and t into a single list for animation purposes
    x_long,y_long,t_long = x[0],y[0],t[0]
    for w in range(1,len(x)):
        x_long = x_long+x[w]
        y_long = y_long+y[w]
        t_long = t_long+t[w]
    
    # rename
    x,y,t = x_long,y_long,t_long
    x,y,t = np.array(x),np.array(y),np.array(t)
    x = max(x)+10-x # ensures that the clock is plotted upright when in 2 dimensions
    
    # order the points by the sum of their x and y coordinates
    sum_xy = []
    for w in range(len(x)):
        sum_xy.append((w,x[w]+y[w]))
    sum_xy = sorted(sum_xy,key=lambda elt: elt[1])
    point_order = [elt[0] for elt in sum_xy]
    x = [x[idx] for idx in point_order]
    y = [y[idx] for idx in point_order]
    t = [t[idx] for idx in point_order]
        
    # plot
    plt.close('all')
    fig_xy = plt.figure()
    fig_xy.text(0.99, 0.96, fname[:len(fname)-4],fontsize=10,color='red',va='baseline',ha='right',multialignment='left')
    xy = fig_xy.add_subplot(111)
    xy.scatter(y,x,color='blue')
    xy.set_xlabel('x',fontsize=20)
    xy.set_ylabel('y',fontsize=20)
    xy.set_xlim(min(y)-10,max(y)+10)
    xy.set_ylim(min(x)-10,max(x)+10)
    plt.axis('equal')
    
    if 'YDU' in fname:
        fig_xy.text(0.25, 0.955, 'HEALTHY',fontsize=15,color='black',va='baseline',ha='right',multialignment='left')
    elif 'CIN' in fname:
        fig_xy.text(0.25, 0.955, 'IMPAIRED',fontsize=15,color='black',va='baseline',ha='right',multialignment='left')
    else:
        print 'not a valid filename'
    
    plt.show()
    
    #fig_xy.savefig(path+'figs_raw/'+fname[:len(fname)-4]+'/whole_'+clock_type+'_clock_angle2d_'+fname[:len(fname)-4]+'.png')
    