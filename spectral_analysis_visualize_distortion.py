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

# path to figures for each individual file
if not os.path.exists(path+'spectral analysis: visualize distortion'):
    os.makedirs(path+'spectral analysis: visualize distortion')

Epeak_vals = [0.8,0.9,1.0]
x_all,y_all = [],[]
N = 250

for Epeak in Epeak_vals:
    dftx,dfty = [0]*250,[0]*250
    dftx[1],dfty[1] = Epeak,Epeak
    dftx[N-1],dfty[N-1] = Epeak,Epeak
    dftx[2],dfty[2] = 1-Epeak,1-Epeak
    dftx[N-2],dfty[N-2] = 1-Epeak,1-Epeak
    
    x = np.fft.ifft(dftx)
    y = np.fft.ifft(dfty)
    y = np.concatenate((y[N/4:],y[:N/4]))
    x_all.append(x)
    y_all.append(y)

for w in range(len(x_all)):
    Epeak_x = Epeak_vals[w]
    x = x_all[w]
    for d in range(w,len(y_all)):
        Epeak_y = Epeak_vals[d]
        y = y_all[d]
    
        # plot for UROP report
        plt.close('all')
        fig_xy = plt.figure()
        xy = fig_xy.add_subplot(111)
        x = [elt-mean(x) for elt in x]
        y = [elt-mean(y) for elt in y]
        xy.plot(x,y,lw=2)
        xy.set_xlim(left=-0.01,right=0.01)
        xy.set_ylim(bottom=-0.01,top=0.01)
    
        # set axis labels
        xy.set_xlabel(r'$x$',fontsize=20)
        xy.set_ylabel(r'$y$',fontsize=20)
        
        j = 1
        for v1 in xy.get_xticklabels():
            v1.set_fontsize(23)
            #if j==1 or j%2==0:
            #    v1.set_visible(False)
            j += 1
        j = 1
        for v2 in xy.get_yticklabels():
            v2.set_fontsize(23)
            #if j==1 or j%2==0:
            #    v2.set_visible(False)
            j += 1
        # save figure
        fig_xy.savefig(path+'spectral analysis: visualize distortion/Epeak_x_'+str(Epeak_x)+'_Epeak_y_'+str(Epeak_y)+'.png')
