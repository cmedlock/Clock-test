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

path = '/Users/cmedlock/Documents/DSP_UROP/simulated data/'

pi = math.pi

# path to normalized data
if not os.path.exists(path+'norm velocity data'):
    os.makedirs(path+'norm velocity data')

# path to figures for each individual file
if not os.path.exists(path+'figs_raw'):
    os.makedirs(path+'figs_raw')

axis_ratios = np.linspace(1,1.2,100)

for ratio in axis_ratios:
    n = np.arange(250)
    x = np.cos(2*np.pi/250*n)
    y = ratio*np.sin(2*np.pi/250*n)
    print 'ratio is now ',ratio
    
    if not os.path.exists(path+'figs_raw/perfect_ellipse_ratio_'+str(ratio)[:5]):
        os.makedirs(path+'figs_raw/perfect_ellipse_ratio_'+str(ratio)[:5])
    
    # compensate for non-constant velocity
    
    N_orig = len(x)
    N_new = N_orig

    # calculate average distance between points
    dists = []
    for w in range(1,len(x)):
        dx,dy = x[w]-x[w-1],y[w]-y[w-1]
        dist = math.sqrt(dx**2+dy**2)
        dists.append(dist)
    dist_total = sum(dists)
    dist_avg = dist_total/float(N_orig)
    #print 'total distance is ',dist_total
    #print 'average distance between points is ',dist_avg

    # now want to get N_orig evenly-spaced points along the curve

    # generate a much longer array with linearly-interpolated 
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
    float_N_new = float(N_new)
    idx = 0
    for k in range(N_new-1):
        dist_sofar = 0
        for j in range(idx,len(x_interp)-1):
            dx,dy = x_interp[j+1]-x_interp[j],y_interp[j+1]-y_interp[j]
            dist_sofar += math.sqrt(dx**2+dy**2)
            if abs(dist_sofar-dist_avg)<dist_avg/100.:
                idx = j+1
                break
        x_eqdist.append(x_interp[idx])
        y_eqdist.append(y_interp[idx])
        
    # subtract mean values so there is no DC term
    x_eqdist = [elt-mean(x_eqdist) for elt in x_eqdist]
    y_eqdist = [elt-mean(y_eqdist) for elt in y_eqdist]
    
    np.savetxt(path+'norm velocity data/perfect_ellipse_x_eqdist_ratio_'+str(ratio)[:5]+'.txt',x_eqdist)
    np.savetxt(path+'norm velocity data/perfect_ellipse_y_eqdist_ratio_'+str(ratio)[:5]+'.txt',y_eqdist)

    # plot for UROP report
    plt.close('all')
    fig_xy = plt.figure()
    xy = fig_xy.add_subplot(111)
    x = [elt-mean(x) for elt in x]
    y = [elt-mean(y) for elt in y]
    xy.plot(x,y,lw=2)
    xy.set_xlim(left=1.1*min(x),right=1.1*max(x))
    xy.set_ylim(bottom=1.1*min(y),top=1.1*max(y))
    
    #frame = plt.gca()
    #frame.axes.get_xaxis().set_ticklabels([])
    #frame.axes.get_yaxis().set_ticklabels([])
    
    # set axis labels
    xy.set_xlabel(r'$x$',fontsize=20)
    xy.set_ylabel(r'$y$',fontsize=20)

    # save figure
    fig_xy.savefig(path+'figs_raw/perfect_ellipse_ratio_'+str(ratio)[:5]+'/xy_raw.png')

    xy.clear()
    xy.plot(x_eqdist,y_eqdist,lw=2)
    xy.set_xlim(left=1.1*min(x_eqdist),right=1.1*max(x_eqdist))
    xy.set_ylim(bottom=1.1*min(y_eqdist),top=1.1*max(y_eqdist))
    
    #frame.axes.get_xaxis().set_ticklabels([])
    #frame.axes.get_yaxis().set_ticklabels([])

    # set axis labels
    xy.set_xlabel(r'$x$',fontsize=20)
    xy.set_ylabel(r'$y$',fontsize=20)

    # save figure
    fig_xy.savefig(path+'figs_raw/perfect_ellipse_ratio_'+str(ratio)[:5]+'/xy.png')
