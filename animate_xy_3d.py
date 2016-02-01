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

# make figs
for fname in dirs:
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

    x,y,t = np.array(x),np.array(y),np.array(t)-t[0]

    # animate
    plt.close('all')
    fig_xy = plt.figure()
    fig_xy.text(0.99, 0.96, fname[:len(fname)-4],fontsize=10,color='red',va='baseline',ha='right',multialignment='left')
    xy = p3.Axes3D(fig_xy)
    xy.set_xlabel('x',fontsize=20)
    xy.set_ylabel('y',fontsize=20)
    xy.set_zlabel('t',fontsize=20)
    xy.set_xlim3d(min(x)-10,max(x)+10)
    xy.set_ylim3d(min(y)-10,max(y)+10)
    xy.set_zlim3d(-5,len(x)+10)

    if 'YDU' in fname:
        fig_xy.text(0.25, 0.955, 'HEALTHY',fontsize=15,color='black',va='baseline',ha='right',multialignment='left')
    elif 'CIN' in fname:
        fig_xy.text(0.25, 0.955, 'IMPAIRED',fontsize=15,color='black',va='baseline',ha='right',multialignment='left')
    else:
        print 'not a valid filename'

    data = [np.array([x,y,np.arange(len(x))])]
    lines = [xy.plot(dat[0, 0:1], dat[1, 0:1], dat[2, 0:1])[0] for dat in data]

    xy_anim = anim.FuncAnimation(fig_xy, update_lines, len(x), fargs=(data, lines), interval=50, blit=False)
    
    plt.draw() # don't delete this
    xy_anim.save(path+'figs_raw/'+fname[:len(fname)-4]+'/'+clock_type+'_clock_3danim_'+fname[:len(fname)-4]+'.mp4',
                      writer=FFwriter,extra_args=['-vcodec','libx264'])
    