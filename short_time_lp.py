# -*- coding: utf-8 -*-
import math
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import numpy.fft
import os
from pylab import *

import clock_test as ct
ct = reload(ct)

path = '/Users/cmedlock/Documents/DSP_UROP/all_data/'
dirs = os.listdir(path)

pi = math.pi

# copy or command clock?
clock_type = 'command'

# path to short time linear prediction figures
if not os.path.exists(path+'short_time_lp'):
    os.makedirs(path+'short_time_lp')
if not os.path.exists(path+'short_time_lp/'+clock_type+'_healthy'):
    os.makedirs(path+'short_time_lp/'+clock_type+'_healthy')
if not os.path.exists(path+'short_time_lp/'+clock_type+'_impaired'):
    os.makedirs(path+'short_time_lp/'+clock_type+'_impaired')

# variance of short time frequency distribution
freq_var_x,freq_var_y = [],[]

for fname in dirs:
    if 'CIN' not in fname and 'YDU' not in fname:
        continue
    #print 'reading file ',fname,'...'
    
    ftype = ''
    if 'YDU' in fname:
        ftype = 'healthy'
    elif 'CIN' in fname:
        ftype = 'impaired'
    else:
        print 'not a valid file name'

    x_eqdist = np.loadtxt(path+'norm_velocity_data/'+fname[:len(fname)-4]+'_x_eqdist_'+clock_type+'.txt')
    y_eqdist = np.loadtxt(path+'norm_velocity_data/'+fname[:len(fname)-4]+'_y_eqdist_'+clock_type+'.txt')

    if len(x_eqdist)==0:
        continue

    # downsample by a factor of 3
    x_downsampled,y_downsampled = [],[]
    for w in range(len(x_eqdist)):
        if w%3==0:
            x_downsampled.append(x_eqdist[w])
            y_downsampled.append(y_eqdist[w])
    #x_eqdist,y_eqdist = np.array(x_downsampled),np.array(y_downsampled)

    # extract frequency from each set of 3 points using the formula
    # w[n] = 2*cos(omega)*w[n-1] - w[n-2]
    x_freq,y_freq = [],[]
    for w in range(len(x_eqdist)):
        frequency_x = (x_eqdist[w]+x_eqdist[w-2])/(2*x_eqdist[w-1])
        frequency_y = (y_eqdist[w]+y_eqdist[w-2])/(2*y_eqdist[w-1])
        x_freq.append(frequency_x)
        y_freq.append(frequency_y)
    freq_var_x.append((ftype,var(x_freq)))
    freq_var_y.append((ftype,var(y_freq)))
    
    # figure settings
    x_axis_fontsize = 15
    y_axis_fontsize = 20

#    # plot
#    plt.close('all')
#    fig = plt.figure()
#    ax = fig.add_subplot(111)
#    ax.scatter(np.arange(len(x_freq)),x_freq,color='blue',marker='o',s=100,alpha=0.5,label='x')
#    ax.scatter(np.arange(len(x_freq)),x_freq,color='red',marker='^',s=100,alpha=0.5,label='y')
#
#    # set axis limits
#    ax.set_ylim(bottom=mean(x_freq)-std(x_freq),top=mean(x_freq)+std(x_freq))
#
#    # set axis labels
#    ax.set_xlabel(r'$n$',fontsize=x_axis_fontsize)
#    ax.set_ylabel(r'cos$(\omega$)',fontsize=y_axis_fontsize)
#
#    # add drawing type (healthy or impaired) and file name
#    fig.text(0.99, 0.96,fname[:len(fname)-4],fontsize=10,color='red',va='baseline',ha='right',multialignment='left')
#
#    if 'YDU' in fname:
#        fig.text(0.32, 0.955, 'HEALTHY ('+clock_type+')',fontsize=15,color='black',va='baseline',ha='right',multialignment='left')
#    elif 'CIN' in fname:
#        fig.text(0.32, 0.955, 'IMPAIRED ('+clock_type+')',fontsize=15,color='black',va='baseline',ha='right',multialignment='left')
#    else:
#        print 'not a valid filename'
#
#    # save figure
#    fig.savefig(path+'figs_raw/'+fname[:len(fname)-4]+'/short_time_lp_'+clock_type+'_'+fname[:len(fname)-4]+'.png')
#    fig.savefig(path+'short_time_lp/'+clock_type+'_'+ftype+'/short_time_lp_'+clock_type+'_'+fname[:len(fname)-4]+'.png')

# compare short time frequency variance for the drawings of healthy vs. impaired patients
freq_var_x_all = [elt[1] for elt in freq_var_x]
mean_x,std_x = mean(freq_var_x_all),std(freq_var_x_all)
print 'mean_x,std_x = ',mean_x,std_x
freq_var_y_all = [elt[1] for elt in freq_var_y]
mean_y,std_y = mean(freq_var_y_all),std(freq_var_y_all)
print 'mean_y,std_y = ',mean_y,std_y
#ct.make_hist([elt[1] for elt in freq_var_x if elt[0]=='healthy'],
#             [elt[1] for elt in freq_var_x if elt[0]=='impaired'],
#             15,mean_x-std_x,mean_x+std_x,'Variance in Short Time Frequency for x[n]','freq_var_x_'+clock_type,path)
#ct.make_hist([elt[1] for elt in freq_var_y if elt[0]=='healthy'],
#             [elt[1] for elt in freq_var_y if elt[0]=='impaired'],
#             15,mean_y-std_y,mean_y+std_y,'Variance in Short Time Frequency for y[n]','freq_var_y_'+clock_type,path)
# in case the histograms don't come out right
np.savetxt(path+'short_time_lp/freq_var_x_healthy_'+clock_type+'.txt',[elt[1] for elt in freq_var_x if elt[0]=='healthy'])
np.savetxt(path+'short_time_lp/freq_var_x_impaired_'+clock_type+'.txt',[elt[1] for elt in freq_var_x if elt[0]=='impaired'])
np.savetxt(path+'short_time_lp/freq_var_y_healthy_'+clock_type+'.txt',[elt[1] for elt in freq_var_y if elt[0]=='healthy'])
np.savetxt(path+'short_time_lp/freq_var_y_impaired_'+clock_type+'.txt',[elt[1] for elt in freq_var_y if elt[0]=='impaired'])
