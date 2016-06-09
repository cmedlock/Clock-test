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
    hist_h,hist_i = fig_hist.add_subplot(211),fig_hist.add_subplot(212)
    n_h,bins_h,patches_h = hist_h.hist(healthy,nbins,range=[range_min,range_max],histtype='bar',color='green',alpha=0.5,label='Healthy')
    n_i,bins_i,patches_i = hist_i.hist(impaired,nbins,range=[range_min,range_max],histtype='bar',color='blue',alpha=0.5,label='Impaired')
    print 'number of entries ignored = ',len(healthy)+len(impaired)-sum(n_h)-sum(n_i),' out of ',len(healthy)+len(impaired)
    hist_h.set_ylim(top=max(max(n_h),max(n_i))*1.2)
    hist_i.set_ylim(top=max(max(n_h),max(n_i))*1.2)

    # set axis label
    x_axis_fontsize = 15
    hist_i.set_xlabel(nameforplot,fontsize=x_axis_fontsize)
    
    j = 1
    for v1 in hist_h.get_yticklabels():
        v1.set_fontsize(23)
        if j==1 or j%2==0:
            v1.set_visible(False)
        j += 1
    j = 1
    for v2 in hist_i.get_yticklabels():
        v2.set_fontsize(23)
        if j==1 or j%2==0:
            v2.set_visible(False)
        j += 1
    j = 1
    for v1 in hist_h.get_xticklabels():
        v1.set_fontsize(17)
        if j%2==0:
            v1.set_visible(False)
        j += 1
    j = 1
    for v2 in hist_i.get_xticklabels():
        v2.set_fontsize(17)
        if j%2==0:
            v2.set_visible(False)
        j += 1
    j = 1

    # save figure
    fig_hist.savefig(path+'/compare_healthy_impaired/compare_'+nameforfile+'.png')
    
    plt.close('all')

def make_hist_xy(healthy_x,impaired_x,healthy_y,impaired_y,nbins,range_min_x,range_max_x,range_min_y,range_max_y,nameforplot,nameforfile,path):
    # input: healthy,impaired (lists): lists of quantity to be compared
    #        name,fname,path: in order to save the figure
    # output: compare_name.png with histogram comparing healthy value to
    #         impaired values

    fig_hist_x = plt.figure()
    hist_h_x,hist_i_x = fig_hist_x.add_subplot(211),fig_hist_x.add_subplot(212)
    n_h,bins_h,patches_h = hist_h_x.hist(healthy_x,nbins,range=[range_min_x,range_max_x],histtype='bar',color='green',alpha=0.5,label='Healthy')
    n_i,bins_i,patches_i = hist_i_x.hist(impaired_x,nbins,range=[range_min_x,range_max_x],histtype='bar',color='blue',alpha=0.5,label='Impaired')
    #print 'number of entries ignored for x[n] = ',len(healthy_x)+len(impaired_x)-sum(n_h)-sum(n_i),' out of ',len(healthy_x)+len(impaired_x)
    print 'healthy patients: mean, std for x[n] = ',mean(healthy_x),std(healthy_x)
    print 'impaired patients: mean, std for x[n] = ',mean(impaired_x),std(impaired_x)
    max_entries_x = max(max(n_h),max(n_i))
    
    fig_hist_y = plt.figure()
    hist_h_y,hist_i_y = fig_hist_y.add_subplot(211),fig_hist_y.add_subplot(212)
    n_h,bins_h,patches_h = hist_h_y.hist(healthy_y,nbins,range=[range_min_y,range_max_y],histtype='bar',color='green',alpha=0.5,label='Healthy')
    n_i,bins_i,patches_i = hist_i_y.hist(impaired_y,nbins,range=[range_min_y,range_max_y],histtype='bar',color='blue',alpha=0.5,label='Impaired')
    #print 'number of entries ignored for y[n] = ',len(healthy_y)+len(impaired_y)-sum(n_h)-sum(n_i),' out of ',len(healthy_y)+len(impaired_y)
    print 'healthy patients: mean, std for y[n] = ',mean(healthy_y),std(healthy_y)
    print 'impaired patients: mean, std for y[n] = ',mean(impaired_y),std(impaired_y)
    max_entries_y = max(max(n_h),max(n_i))
    
    # equalize axis ranges
    max_entries = max(max_entries_x,max_entries_y)
    hist_h_x.set_ylim(top=max_entries*1.2)
    hist_i_x.set_ylim(top=max_entries*1.2)
    hist_h_y.set_ylim(top=max_entries*1.2)
    hist_i_y.set_ylim(top=max_entries*1.2)

    # set axis labels
    x_axis_fontsize = 15
    hist_i_x.set_xlabel(nameforplot,fontsize=x_axis_fontsize)
    hist_i_y.set_xlabel(nameforplot,fontsize=x_axis_fontsize)

    j = 1
    for v1 in hist_h_x.get_yticklabels():
        v1.set_fontsize(23)
        if j==1 or j%2==0:
            v1.set_visible(False)
        j += 1
    j = 1
    for v2 in hist_i_x.get_yticklabels():
        v2.set_fontsize(23)
        if j==1 or j%2==0:
            v2.set_visible(False)
        j += 1
    j = 1
    for v1 in hist_h_x.get_xticklabels():
        v1.set_fontsize(17)
        if j%2==0:
            v1.set_visible(False)
        j += 1
    j = 1
    for v2 in hist_i_x.get_xticklabels():
        v2.set_fontsize(17)
        if j%2==0:
            v2.set_visible(False)
        j += 1
    j = 1

    for v1 in hist_h_y.get_yticklabels():
        v1.set_fontsize(23)
        if j==1 or j%2==0:
            v1.set_visible(False)
        j += 1
    j = 1
    for v2 in hist_i_y.get_yticklabels():
        v2.set_fontsize(23)
        if j==1 or j%2==0:
            v2.set_visible(False)
        j += 1
    j = 1
    for v1 in hist_h_y.get_xticklabels():
        v1.set_fontsize(17)
        if j%2==0:
            v1.set_visible(False)
        j += 1
    j = 1
    for v2 in hist_i_y.get_xticklabels():
        v2.set_fontsize(17)
        if j%2==0:
            v2.set_visible(False)
        j += 1
    j = 1

    # save figure
    fig_hist_x.savefig(path+'/compare_healthy_impaired/compare_'+nameforfile+'_x.png')
    fig_hist_y.savefig(path+'/compare_healthy_impaired/compare_'+nameforfile+'_y.png')
    
    plt.close('all')
