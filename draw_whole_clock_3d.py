# -*- coding: utf-8 -*-
# draw the entire clock, including numbers and clock hands, in 3 dimensions
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
    # order the strokes chronologically
    symbol_start_times = [elt[0] for elt in t]
    symbol_start_times.sort()
    symbol_order = []
    for time in symbol_start_times:
        for w in range(len(t)):
            if t[w][0]==time:
                symbol_order.append(w)
    x = [x[symbol_num] for symbol_num in symbol_order]
    y = [y[symbol_num] for symbol_num in symbol_order]
    t = [t[symbol_num] for symbol_num in symbol_order]
    #for symbol_num in symbol_order:
        #print 'stroke ',stroke_num,': ',t[stroke_num][0],' to ',t[stroke_num][-1]
    
    # get new maxima and minima
    xmin,xmax = min([min(elt) for elt in x]),max([max(elt) for elt in x])
    ymin,ymax = min([min(elt) for elt in y]),max([max(elt) for elt in y])
    tmin,tmax = min([min(elt) for elt in t]),max([max(elt) for elt in t])

    # make each of x, y, and t into a single list for animation purposes
    # insert NaN's between the symbols so that they are not connected together
    x_long,y_long,t_long = x[0]+[np.nan],y[0]+[np.nan],t[0]+[np.nan]
    for w in range(1,len(x)):
        x_long = x_long+x[w]+[np.nan]
        y_long = y_long+y[w]+[np.nan]
        t_long = t_long+t[w]+[np.nan]
    # rename
    x,y,t = np.array(x_long),np.array(y_long),np.array(t_long)

    # plot
    plt.close('all')
    fig_xy = plt.figure()
    fig_xy.text(0.99, 0.96, fname[:len(fname)-4],fontsize=10,color='red',va='baseline',ha='right',multialignment='left')
    xy = fig_xy.add_subplot(111,projection='3d')
    xy.plot(y,xmax+10-x,t,color='blue')
    xy.set_xlabel('y',fontsize=20)
    xy.set_ylabel('x',fontsize=20)
    xy.set_zlabel('t',fontsize=20)
    xy.set_xlim(xmin-10,xmax+10)
    xy.set_ylim(ymin-10,ymax+10)
    xy.set_zlim(tmin-10,tmax+10)
    plt.axis('equal')
    
    if 'YDU' in fname:
        fig_xy.text(0.32, 0.955, 'HEALTHY ('+clock_type+')',fontsize=15,color='black',va='baseline',ha='right',multialignment='left')
    elif 'CIN' in fname:
        fig_xy.text(0.32, 0.955, 'IMPAIRED ('+clock_type+')',fontsize=15,color='black',va='baseline',ha='right',multialignment='left')
    else:
        print 'not a valid filename'
    
    #xy.spines['right'].set_position('zero')
    #xy.spines['top'].set_position('zero')
    #xy.spines['left'].set_position('zero')
    #xy.spines['bottom'].set_position('zero')
    xy.view_init(30,0)
    #plt.show()
    
    fig_xy.savefig(path+'figs_raw/'+fname[:len(fname)-4]+'/whole_'+clock_type+'_clock_3d_'+fname[:len(fname)-4]+'.png')
