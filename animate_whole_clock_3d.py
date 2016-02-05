# -*- coding: utf-8 -*-
# make a video of the entire clock being drawn in 3 dimensions
import math
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as anim
import numpy as np
import os
from pylab import *

def update_lines(num, dataLines, lines):
    for line, data in zip(lines, dataLines):
        # NOTE: there is no .set_data() for 3 dim data...
        line.set_data(data[0:2, :num])
        line.set_3d_properties(data[2, :num])
    return lines

# set-up to save animations
plt.rcParams['animation.ffmpeg_path'] = '/usr/bin/ffmpeg'
FFwriter = anim.FFMpegWriter(fps=1000./13.)

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
    record,found_clock,current_symbol = False,False,''
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
                current_symbol = 'NOISE'
            if 'symbol label' in line and len(x_temp)>0: # length requirement since sometimes there are empty strokes
                #print 'storing previous stroke from ',t_temp[0],' to ',t_temp[-1],' (',(t_temp[0]-1307115563146.0)/10.,' to ',(t_temp[-1]-1307115563146.0)/10.,')'
                #print 'len(t_temp) is now ',len(t_temp)
                #print 'next stroke is ',line
                # store previous stroke and reset lists for the next one
                if current_symbol!='NOISE':
                    x.append(x_temp)
                    y.append(y_temp)
                    t.append(t_temp)
                elif current_symbol=='NOISE':
                    current_symbol = ''
                #print 'len(t) is now ',len(t)
                x_temp,y_temp,t_temp = [],[],[]
                if 'NOISE' in line:
                    current_symbol = 'NOISE'
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
    #for stroke_num in symbol_order:
        #print 'stroke ',stroke_num,': ',t[stroke_num][0],' to ',t[stroke_num][-1]
    
    # make each of x, y, and t into a single list for animation purposes
    x_long,y_long,t_long = x[0],y[0],t[0]
    for w in range(1,len(x)):
        x_long = x_long+x[w]
        y_long = y_long+y[w]
        t_long = t_long+t[w]
    # rename
    x,y,t = x_long,y_long,t_long
    x,y,t = np.array(x),np.array(y),np.array(t)

    # get new maxima and minima
    xmin,xmax = min(x),max(x)
    ymin,ymax = min(y),max(y)
    tmin,tmax = min(t),max(t)
    # plot
    plt.close('all')
    fig_xy = plt.figure()
    fig_xy.text(0.99, 0.96, fname[:len(fname)-4],fontsize=10,color='red',va='baseline',ha='right',multialignment='left')
    #xy = p3.Axes3D(fig_xy)
    xy = fig_xy.add_subplot(111,projection='3d')
    xy.set_xlabel('y',fontsize=20)
    xy.set_ylabel('x',fontsize=20)
    xy.set_zlabel('t',fontsize=20)
    plt.axis('equal')

    if 'YDU' in fname:
        fig_xy.text(0.25, 0.955, 'HEALTHY',fontsize=15,color='black',va='baseline',ha='right',multialignment='left')
    elif 'CIN' in fname:
        fig_xy.text(0.25, 0.955, 'IMPAIRED',fontsize=15,color='black',va='baseline',ha='right',multialignment='left')
    else:
        print 'not a valid filename'

    #data = [np.array([x,y,t])]
    data = [np.array([y,xmax+10-x,t])]
    lines = [xy.plot(dat[0, 0:1], dat[1, 0:1], dat[2, 0:1])[0] for dat in data]

    xy_anim = anim.FuncAnimation(fig_xy, update_lines, len(x), fargs=(data, lines), interval=50, blit=False)
    
    # axis limits need to be set after the animation is drawn (not true for just a 3d figure)
    xy.set_xlim(0,75)
    xy.set_ylim(0,75)
    xy.set_zlim(0,5000)

    xy.view_init(30,0)

    plt.draw() # don't delete this
    xy_anim.save(path+'figs_raw/'+fname[:len(fname)-4]+'/'+clock_type+'_clock_3danim_'+fname[:len(fname)-4]+'.mp4',
                      writer=FFwriter,extra_args=['-vcodec','libx264'])
