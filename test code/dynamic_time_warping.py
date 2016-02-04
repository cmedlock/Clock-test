# calculate DWT between each healthy drawing and each impaired drawing

import math
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
from pylab import *

path = '/Users/cmedlock/Documents/DSP_UROP/subset_data/'
dirs = os.listdir(path)

if not os.path.exists(path+'figs_raw'):
    os.makedirs(path+'figs_raw')

# store signals
healthy_x_copy,healthy_y_copy         = [],[]
impaired_x_copy,impaired_y_copy       = [],[]
healthy_x_command,healthy_y_command   = [],[]
impaired_x_command,impaired_y_command = [],[]

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
        
        # position and time
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
    
    if 'YDU' in fname:
        healthy_x_copy.append(x_copy)
        healthy_y_copy.append(y_copy)
        healthy_x_command.append(x_command)
        healthy_y_command.append(y_command)
    elif 'CIN' in fname:
        impaired_x_copy.append(x_copy)
        impaired_y_copy.append(y_copy)
        impaired_x_command.append(x_command)
        impaired_y_command.append(y_command)

for w in range(len(healthy_x_copy)):
    for d in range(len(impaired_x_copy)):
        #Q,n = healthy_x_copy[w],len(healthy_x_copy[w])
        #C,m = impaired_x_copy[d],len(impaired_x_copy[d])
        Q,n = range(4),4
        C,m = range(5),5
        print 'n = ',n,' and m = ',m
        # DTW matrix
        D = np.zeros((n,m))
        # fill it in and find warping path
        W = [(0,0)]
        for v in range(n):
            for x in range(m):
                cost = (Q[v]-C[x])**2
                print '   v = ',v,', x = ',x,', and cost = ',cost
                D[v][x] = cost + min(D[v-1][x],D[v][x-1],D[v-1,x-1])
        print 'DTW is ',D[n-1][m-1]
