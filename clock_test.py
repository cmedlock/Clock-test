import math
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as anim
import numpy as np
import numpy.fft
import os
from pylab import *

def deriv_doublederiv(x,y,t):
    # input: x[t],y[t],t
    # output: vx[t],vy[t],ax[t],ay[t]
    
    if len(x)!=len(y):
        raise NameError('ERROR: x and y lists do not have the same length')

    vx,vy = [],[]
    for d in range(1,len(x)):
        v_x = (x[d]-x[d-1])/(t[d]-t[d-1])
        v_y = (y[d]-y[d-1])/(t[d]-t[d-1])
        vx.append(v_x)
        vy.append(v_y)
    ax,ay = [],[]
    for b in range(1,len(vx)):
        a_x = (vx[b]-vx[b-1])/(t[b]-t[b-1])
        a_y = (vy[b]-vy[b-1])/(t[b]-t[b-1])
        ax.append(a_x)
        ay.append(a_y)
    return vx,vy,ax,ay

def get_bins(a,nbins):
    # input: a (list of tuples): [('healthy',...),('impaired',...)...]
    #       want to find appropriate bin edges so that can make a
    #       histogram of the healthy vs. the impaired with equal bin widths
    # output: binedges (list of tuples with the bin edges)
    healthy = [elt[1] for elt in a if elt[0]=='healthy']
    impaired = [elt[1] for elt in a if elt[0]=='impaired']
    binedges = np.histogram(np.hstack((healthy,impaired)),bins=10)[1]
    return binedges
    
def make_hist(healthy,impaired,nbins,range_min,range_max,nameforplot,nameforfile,path):
    # input: healthy,impaired (lists): lists of quantity to be compared
    #        name,fname,path: in order to save the figure
    # output: compare_name.png with histogram comparing healthy value to
    #         impaired values

    fig_hist = plt.figure()
    hist = fig_hist.add_subplot(111)
    n,bins,patches = hist.hist([healthy,impaired],nbins,range=[range_min,range_max],histtype='bar',label=['Healthy','Impaired'])
    n1,n2 = n[0],n[1]
    print 'number of entries ignored = ',len(healthy)+len(impaired)-sum(n1)-sum(n2),' out of ',len(healthy)+len(impaired)
    hist.set_ylim(top=max(max(n1),max(n2))*1.2)
    hist.legend(loc='best',frameon=False)
    
    # set axis label
    x_axis_fontsize = 15
    hist.set_xlabel(nameforplot,fontsize=x_axis_fontsize)

    # save figure
    fig_hist.savefig(path+'/compare_healthy_impaired/compare_'+nameforfile+'.png')
    
    plt.close('all')
    