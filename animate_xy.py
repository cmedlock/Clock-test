import math
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as anim
import numpy as np
import os
from pylab import *

def init_anim():
    line.set_data([],[])
    return line,
    
def animate(nitr,a,b):
    line.set_data(a[:nitr+1],b[:nitr+1])
    return line,

# set-up to save animations
plt.rcParams['animation.ffmpeg_path'] = '/usr/bin/ffmpeg'
FFwriter = anim.FFMpegWriter(fps=1000./13.)

# path to data
path = '/Users/cmedlock/Documents/DSP UROP/all_data/'
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
    clock_type = ''
    x_copy,y_copy,t_copy = [],[],[]
    x_command,y_command,t_command = [],[],[]

    # read in data
    record,found_clock = False,False
    for w in range(len(data)):
        line = data[w]
        # found copy clock?
        if found_clock==False:
            if 'COPY' in line:
                clock_type = 'COPY'
                found_clock = True
            elif 'COMMAND' in line:
                clock_type = 'COMMAND'
                found_clock = True
            continue
        # start recording?
        elif found_clock==True and 'CLOCKFACE' in line:
            record = True
            continue
        # stop recording?
        elif record==True:
            if 'symbol label' in line and clock_type=='COPY' and len(x_copy)>0:
                found_clock = False
                record = False
                continue
            elif 'symbol label' in line and clock_type=='COMMAND' and len(x_command)>0:
                found_clock = False
                record = False
                continue
            elif 'point' not in line:
                continue
        # other?
        else:
            continue
        line = line.split('"')

        # get position and time
        xcoord = double(line[3])
        ycoord = double(line[1])
        timestamp = double(line[7])
        if clock_type=='COPY':
            x_copy.append(xcoord)
            y_copy.append(ycoord)
            t_copy.append(timestamp)
        elif clock_type=='COMMAND':
            x_command.append(xcoord)
            y_command.append(ycoord)
            t_command.append(timestamp)
        else:
            print 'not a valid clock type'
    
    f.close()

    x_copy,y_copy,t_copy = np.array(x_copy),np.array(y_copy),np.array(t_copy)
    x_command,y_command,t_command = np.array(x_command),np.array(y_command),np.array(t_command)

    # copy clock
    fig_xy_copy = plt.figure()
    fig_xy_copy.text(0.99, 0.96, fname[:len(fname)-4],fontsize=10,color='red',va='baseline',ha='right',multialignment='left')
    xy_copy = fig_xy_copy.add_subplot(111,aspect=1.0)
    xy_copy.set_xlabel('x',fontsize=20)
    xy_copy.set_ylabel('y',fontsize=20)
    xy_copy.set_xlim(min(x_copy)-10,max(x_copy)+10)
    xy_copy.set_ylim(min(y_copy)-10,max(y_copy)+10)

    line, = xy_copy.plot([],[],lw=2)
    xy_copy_anim = anim.FuncAnimation(fig_xy_copy,animate,frames=len(x_copy),
                                      init_func=init_anim,fargs=(x_copy,y_copy),
                                      interval=13,blit=True,repeat=False)                     
    plt.draw() # don't delete this
    xy_copy_anim.save(path+'figs_raw/'+fname[:len(fname)-4]+'/copy_clock_anim_'+fname[:len(fname)-4]+'.mp4',
                      writer=FFwriter,extra_args=['-vcodec','libx264'])
    
    # command clock
    fig_xy_command = plt.figure()
    fig_xy_command.text(0.99, 0.96, fname[:len(fname)-4],fontsize=10,color='red',va='baseline',ha='right',multialignment='left')
    xy_command = fig_xy_command.add_subplot(111,aspect=1.0)
    xy_command.set_xlabel('x',fontsize=20)
    xy_command.set_ylabel('y',fontsize=20)
    xy_command.set_xlim(min(x_command)-10,max(x_command)+10)
    xy_command.set_ylim(min(y_command)-10,max(y_command)+10)
    line, = xy_command.plot([],[],lw=2)
    xy_command_anim = anim.FuncAnimation(fig_xy_command,animate,frames=len(x_command),
                                      init_func=init_anim,fargs=(x_command,y_command),
                                      interval=13,blit=True,repeat=False)
    plt.draw() # don't delete this
    xy_command_anim.save(path+'figs_raw/'+fname[:len(fname)-4]+'/command_clock_anim_'+fname[:len(fname)-4]+'.mp4',
                         writer=FFwriter,extra_args=['-vcodec','libx264'])
 
    plt.close('all')
