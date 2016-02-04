# -*- coding: utf-8 -*-
# draw the entire clock, including numbers and clock hands
import math
import matplotlib
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as anim
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
    clock_type = 'COPY'
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
            if 'symbol label' in line and len(x_temp)>0: # length requirement since sometimes there are empty strokes
                #print 'storing previous stroke from ',t_temp[0],' to ',t_temp[-1],' (',(t_temp[0]-1307115563146.0)/10.,' to ',(t_temp[-1]-1307115563146.0)/10.,')'
                #print 'len(t_temp) is now ',len(t_temp)
                #print 'next stroke is ',line
                # store previous stroke and reset lists for the next one
                if is_noise==False:
                    x.append(x_temp)
                    y.append(y_temp)
                    t.append(t_temp)
                elif is_noise==True:
                    is_noise = False
                print 'len(t) is now ',len(t)
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
    x = [np.array(elt)-xmin for elt in x]
    y = [np.array(elt)-ymin for elt in y]
    t = [(np.array(elt)-tmin)/10. for elt in t] # the time array should have reasonable values
    # order the strokes chronologically
    stroke_start_times = [elt[0] for elt in t]
    stroke_start_times.sort()
    stroke_order = []
    for time in stroke_start_times:
        for w in range(len(t)):
            if t[w][0]==time:
                stroke_order.append(w)
    for stroke_num in stroke_order:
        print 'stroke ',stroke_num,': ',t[stroke_num][0],' to ',t[stroke_num][-1]

    xmax = max([max(elt) for elt in x])
    # plot
    plt.close('all')
    fig_xy = plt.figure()
    fig_xy.text(0.99, 0.96, fname[:len(fname)-4],fontsize=10,color='red',va='baseline',ha='right',multialignment='left')
    xy = fig_xy.add_subplot(111)
    for w in range(len(x)):
        xy.plot(y[w],xmax+10-x[w],color='blue')
    xy.set_xlabel('x',fontsize=20)
    xy.set_ylabel('y',fontsize=20)
    xy.set_xlim(min(x[0])-10,max(x[0])+10)
    xy.set_ylim(min(y[0])-10,max(y[0])+10)
    plt.axis('equal')
    
    if 'YDU' in fname:
        fig_xy.text(0.25, 0.955, 'HEALTHY',fontsize=15,color='black',va='baseline',ha='right',multialignment='left')
    elif 'CIN' in fname:
        fig_xy.text(0.25, 0.955, 'IMPAIRED',fontsize=15,color='black',va='baseline',ha='right',multialignment='left')
    else:
        print 'not a valid filename'
    
    fig_xy.savefig(path+'figs_raw/'+fname[:len(fname)-4]+'/whole_'+clock_type+'_clock_'+fname[:len(fname)-4]+'.png')
