# -*- coding: utf-8 -*-
# make a video of the clock being drawn, from an angle different than the
# most straightforward one (looking straight down the time axis at the xy
# plane
import math
import matplotlib
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
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

# copy or command clock?
clock_type = 'COMMAND'

# make figs
for fname in dirs[4:5]:
    if 'Scored' not in fname:
        continue
    print 'reading file ',fname,'...'
    f = open(path+fname)
    data = f.readlines()

    if not os.path.exists(path+'figs_raw/'+fname[:len(fname)-4]):
        os.makedirs(path+'figs_raw/'+fname[:len(fname)-4])

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

    # change x so that if the whole clock is drawn,
    # it is oriented correctly
    xmax = max([max(elt) for elt in x])
    for w in range(len(x)):
        x[w] = [xmax+10-elt for elt in x[w]]
    xmax = max([max(elt) for elt in x])
    xmin = min([min(elt) for elt in x])
    ymax = max([max(elt) for elt in y])
    ymin = min([min(elt) for elt in y])
    tmin = min([min(elt) for elt in t])
    # also normalize the size of the word such that all coordinates
    # are linearly positioned between 0 and 127
    for w in range(len(x)):
        x[w] = [127*(elt-xmin)/(xmax-xmin) for elt in x[w]]
        y[w] = [127*(elt-ymin)/(ymax-ymin) for elt in y[w]]
        t[w] = [(elt-tmin)/1000. for elt in t[w]] # the time array should have reasonable values
    
    # label each point with a number that represents what symbol it is a part of
    # this is necessary for animating the clock at a different angle, since the symbols
    # won't be drawn one after another anymore
    for w in range(len(x)):
        for d in range(len(x[w])):
            x[w][d] = (w,x[w][d])
            y[w][d] = (w,y[w][d])
            t[w][d] = (w,t[w][d])
           
    # make each of x, y, and t into a single list for animation purposes
    # insert NaN's between the symbols so that they are not connected together
    x_long,y_long,t_long = x[0],y[0],t[0]
    for w in range(1,len(x)):
        x_long = x_long+x[w]
        y_long = y_long+y[w]
        t_long = t_long+t[w]
    # rename
    x,y,t = x_long,y_long,t_long

    # order the points by the sum of their x and y coordinates
    sum_xy = []
    for w in range(len(x)):
        sum_xy.append((w,x[w][1]+y[w][1]))
        #sum_xy.append((w,x[w][1]+y[w][1]+5*t[w][1]))
    sum_xy = sorted(sum_xy,key=lambda elt: elt[1])
    point_order = [elt[0] for elt in sum_xy]
    x = [x[idx][1] for idx in point_order][::-1]
    y = [y[idx][1] for idx in point_order][::-1]
    t = [t[idx][1] for idx in point_order][::-1]

    # plot
    plt.close('all')
    fig_xy = plt.figure()
    fig_xy.text(0.99, 0.96, fname[:len(fname)-4],fontsize=10,color='red',va='baseline',ha='right',multialignment='left')
    xy = fig_xy.add_subplot(111,projection='3d')
    #xy.plot(y,x,t,color='blue')
    xy.set_xlabel('x',fontsize=20)
    xy.set_ylabel('y',fontsize=20)
    xy.set_zlabel('t',fontsize=20)
    
    if 'YDU' in fname:
        fig_xy.text(0.32, 0.955, 'HEALTHY ('+clock_type,+')',fontsize=15,color='black',va='baseline',ha='right',multialignment='left')
    elif 'CIN' in fname:
        fig_xy.text(0.32, 0.955, 'IMPAIRED ('+clock_type+')',fontsize=15,color='black',va='baseline',ha='right',multialignment='left')
    else:
        print 'not a valid filename'

    data = [np.array([y,x,t])] # connect all symbols
    lines = [xy.plot(dat[0, 0:1], dat[1, 0:1], dat[2, 0:1],'o',markersize=2)[0] for dat in data]

    xy_anim = anim.FuncAnimation(fig_xy, update_lines, len(x), fargs=(data, lines), interval=13, blit=False) # connect all symbols

    xy.view_init(0,-45) # check the regular animation (projection onto the xy-plane)
    
    # axis limits need to be set after the animation is drawn (not true for a 3d figure)
    xy.set_xlim(-10,140)
    xy.set_ylim(-10,140)
    xy.set_zlim(-10,45)

    plt.draw() # don't delete this
    #xy_anim.save(path+'figs_raw/'+fname[:len(fname)-4]+'/'+clock_type+'_whole_clock_anim_angle2d_'+fname[:len(fname)-4]+'.mp4',
    #             writer=FFwriter,extra_args=['-vcodec','libx264'])

    plt.show()
    
    #fig_xy.savefig(path+'figs_raw/'+fname[:len(fname)-4]+'/whole_'+clock_type+'_clock_angle2d_'+fname[:len(fname)-4]+'.png')
