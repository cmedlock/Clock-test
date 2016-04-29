# plot DFS coefficients of each x vs. t and y vs. t signal

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

# copy or command clock?
clock_type = 'command'

# path to DFS coefficient figures
if not os.path.exists(path+'DFS_coefficients'):
    os.makedirs(path+'DFS_coefficients')
if not os.path.exists(path+'DFS_coefficients/'+clock_type+'_healthy'):
    os.makedirs(path+'DFS_coefficients/'+clock_type+'_healthy')
if not os.path.exists(path+'DFS_coefficients/'+clock_type+'_impaired'):
    os.makedirs(path+'DFS_coefficients/'+clock_type+'_impaired')
if not os.path.exists(path+'compare_healthy_impaired'):
    os.makedirs(path+'compare_healthy_impaired')

# fraction of energy contained in fundamental frequency
Epeak_x,Epeak_y = [],[]

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
    
    # DFS coefficients
    dft_size = len(x_eqdist)
    dftx,dfty = np.fft.fft(x_eqdist,n=dft_size),np.fft.fft(y_eqdist,n=dft_size)
    k = np.arange(dft_size)

    # k_near_pi is the k value for which w_k = 2*pi*k/N is closest to,
    # but not larger than, pi
    k_near_pi = 0
    if dft_size%2==0:
        k_near_pi = dft_size/2+1
    else:
        k_near_pi = math.ceil(dft_size/2)

    # percent energy in peak
    Ex,Ey = np.abs(dftx)**2,np.abs(dfty)**2
    Ex_total,Ey_total = sum(Ex),sum(Ey)
    Ex_peak,Ey_peak = 2*Ex[1]/Ex_total,2*Ey[1]/Ey_total
    Epeak_x.append((ftype,Ex_peak))
    Epeak_y.append((ftype,Ey_peak))

    # only use the positive frequencies
    posfreq = k[:k_near_pi]
    dftx_posfreq,dfty_posfreq = dftx[:k_near_pi],dfty[:k_near_pi]
    #freq = 2*np.pi/dft_size*k
    #freq = np.fft.fftfreq(n=1024,d=1/(2*np.pi)) # NB: centered around 0
    
#    # plot
#    fig_xt,fig_yt = plt.figure(),plt.figure()
#    fig_xt.subplots_adjust(hspace=0.3)
#    fig_yt.subplots_adjust(hspace=0.3)
#    
#    dftxk,dftyk = fig_xt.add_subplot(211),fig_yt.add_subplot(211)    
#    linecolor_dft = 'red'
#    dftxk.plot(posfreq,np.abs(dftx_posfreq),linecolor_dft)
#    dftyk.plot(posfreq,np.abs(dfty_posfreq),linecolor_dft)
#    dftxk.set_ylim(bottom=min(np.abs(dftx_posfreq[1:]))/10,top=max(np.abs(dftx))*10)
#    dftyk.set_ylim(bottom=min(np.abs(dfty_posfreq[1:]))/10,top=max(np.abs(dfty))*10)
#    dftxk.set_yscale('log')
#    dftyk.set_yscale('log')
#
#    dftxk_zoom,dftyk_zoom = fig_xt.add_subplot(212),fig_yt.add_subplot(212)
#    dftxk_zoom.stem(posfreq[:10],np.abs(dftx_posfreq[:10]),linecolor_dft,markerfmt='ro')
#    dftyk_zoom.stem(posfreq[:10],np.abs(dfty_posfreq[:10]),linecolor_dft,markerfmt='ro')
#    dftxk_zoom.set_ylim(bottom=min(np.abs(dftx_posfreq[1:10]))/10,top=max(np.abs(dftx))*10)
#    dftyk_zoom.set_ylim(bottom=min(np.abs(dfty_posfreq[1:10]))/10,top=max(np.abs(dfty))*10)
#    dftxk_zoom.set_yscale('log')
#    dftyk_zoom.set_yscale('log')
#    
#    # set axis limits
#    dftxk.set_xlim(right=posfreq[-1])
#    dftyk.set_xlim(right=posfreq[-1])
#    dftxk_zoom.set_xlim(left=0,right=posfreq[10])
#    dftyk_zoom.set_xlim(left=0,right=posfreq[10])
#    
#    # set axis labels
#    x_axis_fontsize = 20
#    dftxk.set_xlabel(r'$k$',fontsize=x_axis_fontsize)
#    dftyk.set_xlabel(r'$k$',fontsize=x_axis_fontsize)
#    dftxk_zoom.set_xlabel(r'$k$',fontsize=x_axis_fontsize)
#    dftyk_zoom.set_xlabel(r'$k$',fontsize=x_axis_fontsize)
#
#    y_axis_fontsize = 25
#    linecolor_xy = 'blue'
#        
#    dftxk.set_ylabel(r'$|X[k]|$',color=linecolor_dft,fontsize=y_axis_fontsize)
#    dftyk.set_ylabel(r'$|Y[k]|$',color=linecolor_dft,fontsize=y_axis_fontsize)
#    for v1 in dftxk.get_yticklabels():
#        v1.set_color(linecolor_dft)
#    for v2 in dftyk.get_yticklabels():
#        v2.set_color(linecolor_dft)
#    dftxk_zoom.set_ylabel(r'$|X[k]|$',color=linecolor_dft,fontsize=y_axis_fontsize)
#    dftyk_zoom.set_ylabel(r'$|Y[k]|$',color=linecolor_dft,fontsize=y_axis_fontsize)
#    for v1 in dftxk_zoom.get_yticklabels():
#        v1.set_color(linecolor_dft)
#    for v2 in dftyk_zoom.get_yticklabels():
#        v2.set_color(linecolor_dft)
#    
#    # add drawing type (healthy or impaired) and file name
#    fig_xt.text(0.99, 0.96,fname[:len(fname)-4],fontsize=10,color='red',va='baseline',ha='right',multialignment='left')
#    fig_yt.text(0.99, 0.96,fname[:len(fname)-4],fontsize=10,color='red',va='baseline',ha='right',multialignment='left')
#
#    if 'YDU' in fname:
#        fig_xt.text(0.25, 0.955, 'HEALTHY',fontsize=15,color='black',va='baseline',ha='right',multialignment='left')
#        fig_yt.text(0.25, 0.955, 'HEALTHY',fontsize=15,color='black',va='baseline',ha='right',multialignment='left')
#    elif 'CIN' in fname:
#        fig_xt.text(0.25, 0.955, 'IMPAIRED',fontsize=15,color='black',va='baseline',ha='right',multialignment='left')
#        fig_yt.text(0.25, 0.955, 'IMPAIRED',fontsize=15,color='black',va='baseline',ha='right',multialignment='left')
#    else:
#        print 'not a valid filename'
#        
#    # save figures
#    fig_xt.savefig(path+'figs_raw/'+fname[:len(fname)-4]+'/dftx_'+clock_type+'_'+fname[:len(fname)-4]+'.png')
#    fig_xt.savefig(path+'DFS_coefficients/'+clock_type+'_'+ftype+'/dftx_'+clock_type+'_'+fname[:len(fname)-4]+'.png')
#    fig_yt.savefig(path+'figs_raw/'+fname[:len(fname)-4]+'/dfty_'+clock_type+'_'+fname[:len(fname)-4]+'.png')
#    fig_yt.savefig(path+'DFS_coefficients/'+clock_type+'_'+ftype+'/dfty_'+clock_type+'_'+fname[:len(fname)-4]+'.png')

    plt.close('all')

# compare energy in fundamental frequency for the drawings of healthy vs. impaired patients
Epeak_x_all = [elt[1] for elt in Epeak_x]
mean_x,std_x = mean(Epeak_x_all),std(Epeak_x_all)
Epeak_y_all = [elt[1] for elt in Epeak_y]
mean_y,std_y = mean(Epeak_y_all),std(Epeak_y_all)
ct.make_hist([elt[1] for elt in Epeak_x if elt[0]=='healthy'],
             [elt[1] for elt in Epeak_x if elt[0]=='impaired'],
             15,mean_x-std_x/5,1.001,'Fraction of Energy in Fundamental Freq. for x[n]','Epeak_x_'+clock_type,path)
ct.make_hist([elt[1] for elt in Epeak_y if elt[0]=='healthy'],
             [elt[1] for elt in Epeak_y if elt[0]=='impaired'],
             15,mean_y-std_y/5,1.001,'Fraction of Energy in Fundamental Freq. for y[n]','Epeak_y_'+clock_type,path)
# in case the histograms don't come out right
np.savetxt(path+'DFS_coefficients/Epeak_x_healthy_'+clock_type+'.txt',[elt[1] for elt in Epeak_x if elt[0]=='healthy'])
np.savetxt(path+'DFS_coefficients/Epeak_x_impaired_'+clock_type+'.txt',[elt[1] for elt in Epeak_x if elt[0]=='impaired'])
np.savetxt(path+'DFS_coefficients/Epeak_y_healthy_'+clock_type+'.txt',[elt[1] for elt in Epeak_y if elt[0]=='healthy'])
np.savetxt(path+'DFS_coefficients/Epeak_y_impaired_'+clock_type+'.txt',[elt[1] for elt in Epeak_y if elt[0]=='impaired'])
