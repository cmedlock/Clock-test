# 1.normalize for non-constant velocity
# 2.estimate ideal underlying sinusoids
# 3.compute correlation and mean squared difference 
# between ideal sinusoids and actual drawing
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

if not os.path.exists(path+'figs_raw'):
    os.makedirs(path+'figs_raw')

# save interesting quantities
corr_x_copy,corr_y_copy = [],[]
corr_x_command,corr_y_command = [],[]
Ediff_x_copy,Ediff_y_copy = [],[]
Ediff_x_command,Ediff_y_command = [],[]

for fname in dirs[:3]:
    if 'Scored' not in fname:
        continue
    print 'reading file ',fname,'...'
    
    ftype = ''
    if 'YDU' in fname:
        ftype = 'healthy'
    elif 'CIN' in fname:
        ftype = 'impaired'
    else:
        print 'not a valid file name'

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

    x_copy,y_copy,t_copy = np.array(x_copy),np.array(y_copy),np.array(t_copy)
    x_command,y_command,t_command = np.array(x_command),np.array(y_command),np.array(t_command)
    
    # compensate for non-constant velocity
    N_orig_copy,N_orig_command = len(x_copy),len(x_command)

    # calculate average distance between points
    # copy clocks
    dists_copy = []
    for w in range(1,len(x_copy)):
        dx,dy = x_copy[w]-x_copy[w-1],y_copy[w]-y_copy[w-1]
        dist = math.sqrt(dx**2+dy**2)
        dists_copy.append(dist)
    dist_avg_copy = mean(dists_copy)
    print 'average distance between points is ',dist_avg_copy
    print 'total distance is ',sum(dists_copy)
    # command clocks
    dists_command = []
    for w in range(1,len(x_command)):
        dx,dy = x_command[w]-x_command[w-1],y_command[w]-y_command[w-1]
        dist = math.sqrt(dx**2+dy**2)
        dists_command.append(dist)
    dist_avg_command = mean(dists_command)

    # now want to get N_orig evenly-spaced points along the curve

    # generate a much longer array with 199 linearly-interpolated 
    # points between the actual data points
    # copy clocks
    x_interp_copy,y_interp_copy = [],[]
    for w in range(len(x_copy)-1):
        x_interp_copy.append(x_copy[w])
        y_interp_copy.append(y_copy[w])
        dx,dy = x_copy[w+1]-x_copy[w],y_copy[w+1]-y_copy[w]
        dist = math.sqrt(dx**2+dy**2)
        for r in range(1,200):
            x_new = x_copy[w]+r*dx/200
            y_new = y_copy[w]+r*dy/200
            x_interp_copy.append(x_new)
            y_interp_copy.append(y_new)
    x_interp_copy.append(x_copy[-1])
    y_interp_copy.append(y_copy[-1])
    #print '\naverage distance between interpolated points is ',dist_avg_interp
    #print 'total distance is now ',sum(dists_interp)
    # command clocks
    x_interp_command,y_interp_command = [],[]
    for w in range(len(x_command)-1):
        x_interp_command.append(x_command[w])
        y_interp_command.append(y_command[w])
        dx,dy = x_command[w+1]-x_command[w],y_command[w+1]-y_command[w]
        dist = math.sqrt(dx**2+dy**2)
        for r in range(1,200):
            x_new = x_command[w]+r*dx/200
            y_new = y_command[w]+r*dy/200
            x_interp_command.append(x_new)
            y_interp_command.append(y_new)
    x_interp_command.append(x_command[-1])
    y_interp_command.append(y_command[-1])
    #print '\naverage distance between interpolated points is ',dist_avg_interp
    #print 'total distance is now ',sum(dists_interp)

    # start from the first point and find the ones that are 
    # approximately a distance dist_avg from each other
    # copy clocks
    x_eqdist_copy,y_eqdist_copy = [x_interp_copy[0]],[y_interp_copy[0]]
    idx = 0
    for k in range(len(x_copy)-1):
        dist_total = 0
        for j in range(idx,len(x_interp_copy)-1):
            dx,dy = x_interp_copy[j+1]-x_interp_copy[j],y_interp_copy[j+1]-y_interp_copy[j]
            dist_total += math.sqrt(dx**2+dy**2)
            if abs(dist_total-dist_avg_copy)<0.01:
                idx = j+1
	        break
        x_eqdist_copy.append(x_interp_copy[idx])
        y_eqdist_copy.append(y_interp_copy[idx])
    x_eqdist_copy,y_eqdist_copy = np.array(x_eqdist_copy),np.array(y_eqdist_copy)
    # command clocks
    x_eqdist_command,y_eqdist_command = [x_interp_command[0]],[y_interp_command[0]]
    idx = 0
    for k in range(len(x_command)-1):
        dist_total = 0
        for j in range(idx,len(x_interp_command)-1):
            dx,dy = x_interp_command[j+1]-x_interp_command[j],y_interp_command[j+1]-y_interp_command[j]
            dist_total += math.sqrt(dx**2+dy**2)
            if abs(dist_total-dist_avg_command)<0.01:
                idx = j+1
	        break
        x_eqdist_command.append(x_interp_command[idx])
        y_eqdist_command.append(y_interp_command[idx])
    x_eqdist_command,y_eqdist_command = np.array(x_eqdist_command),np.array(y_eqdist_command)

    # now want to estimate the frequency of the underlying sinusoid
    # subtract mean values (zm = zero mean)
    x_eqdist_zm_copy,y_eqdist_zm_copy = x_eqdist_copy-mean(x_eqdist_copy),y_eqdist_copy-mean(y_eqdist_copy)
    x_eqdist_zm_command,y_eqdist_zm_command = x_eqdist_command-mean(x_eqdist_command),y_eqdist_command-mean(y_eqdist_command)
    
    # DFS coefficients
    dft_size_copy,dft_size_command = N_orig_copy,N_orig_command
    k_copy,k_command = np.arange(dft_size_copy),np.arange(dft_size_command)
    dftx_copy,dfty_copy = np.fft.fft(x_eqdist_zm_copy,n=dft_size_copy),np.fft.fft(y_eqdist_zm_copy,n=dft_size_copy)
    dftx_command,dfty_command = np.fft.fft(x_eqdist_zm_command,n=dft_size_command),np.fft.fft(y_eqdist_zm_command,n=dft_size_command)

    # k_near_pi is the k value for which w_k = 2*pi*k/N is closest to,
    # but not larger than, pi
    k_near_pi_copy,k_near_pi_command = 0,0
    if dft_size_copy%2==0:
        k_near_pi_copy = dft_size_copy/2+1
    else:
        k_near_pi_copy = math.ceil(dft_size_copy/2)
    if dft_size_command%2==0:
        k_near_pi_command = dft_size_command/2+1
    else:
        k_near_pi_command = math.ceil(dft_size_command/2)

    # only use the positive frequencies
    posfreq_copy = k_copy[:k_near_pi_copy]
    posfreq_command = k_command[:k_near_pi_command]
    dftx_posfreq_copy,dfty_posfreq_copy = dftx_copy[:k_near_pi_copy],dfty_copy[:k_near_pi_copy]
    dftx_posfreq_command,dfty_posfreq_command = dftx_command[:k_near_pi_command],dfty_command[:k_near_pi_command]
    #freq = 2*np.pi/dft_size*k
    #freq = np.fft.fftfreq(n=1024,d=1/(2*np.pi)) # NB: centered around 0
    
    # take the frequency of the largest DFS coefficient to be
    # the approximate frequency of the underlying sinusoid
    # copy clocks
    abs_dftx_posfreq_copy,abs_dfty_posfreq_copy = list(np.abs(dftx_posfreq_copy)),list(np.abs(dfty_posfreq_copy))
    k_true_x_copy = abs_dftx_posfreq_copy.index(max(abs_dftx_posfreq_copy))
    k_true_y_copy = abs_dfty_posfreq_copy.index(max(abs_dfty_posfreq_copy))
    w_true_x_copy = 2*pi*k_true_x_copy/dft_size_copy
    w_true_y_copy = 2*pi*k_true_y_copy/dft_size_copy
    #print '\nw_true_x_copy = ',w_true_x_copy,' and w_true_y_copy = ',w_true_y_copy
    #print '\n',np.angle(dftx[1]),np.angle(dftx[-1]),' ---> ',np.angle(dftx[1])-np.angle(dftx[-1])
    #print np.angle(dfty[1]),np.angle(dfty[-1]),' ---> ',np.angle(dfty[1])-np.angle(dfty[-1])
    # command clocks
    abs_dftx_posfreq_command,abs_dfty_posfreq_command = list(np.abs(dftx_posfreq_command)),list(np.abs(dfty_posfreq_command))
    k_true_x_command = abs_dftx_posfreq_command.index(max(abs_dftx_posfreq_command))
    k_true_y_command = abs_dfty_posfreq_command.index(max(abs_dfty_posfreq_command))
    w_true_x_command = 2*pi*k_true_x_command/dft_size_command
    w_true_y_command = 2*pi*k_true_y_command/dft_size_command

    # calculate ideal underlying sinusoid
    n_copy,n_command = np.arange(len(x_eqdist_copy)),np.arange(len(x_eqdist_command))
    x_true_copy,y_true_copy = np.cos(w_true_x_copy*n_copy),np.sin(w_true_y_copy*n_copy)
    x_true_command,y_true_command = np.cos(w_true_x_command*n_command),np.sin(w_true_y_command*n_command)
    # use maximum correlation to determine phase and range of x or y to determine amplitude
    # copy clocks
    phase_x_copy,max_corr_x_copy = 0,0
    phase_y_copy,max_corr_y_copy = 0,0
    for w in range(len(x_true_copy)):
        x_shifted = np.concatenate((x_true_copy[w:],x_true_copy[:w]))
        y_shifted = np.concatenate((y_true_copy[w:],y_true_copy[:w]))
        corr_x = np.dot(x_eqdist_copy,x_shifted)
        corr_y = np.dot(y_eqdist_copy,y_shifted)
        if corr_x>max_corr_x_copy:
            max_corr_x_copy = corr_x
            phase_x_copy = w
        if corr_y>max_corr_y_copy:
            max_corr_y_copy = corr_y
            phase_y_copy = w
    #print 'phase for x_true[n] is ',phase_x,' and max correlation is ',max_corr_x
    #print 'phase for y_true[n] is ',phase_y,' and max correlation is ',max_corr_y
    x_true_copy = np.concatenate((x_true_copy[phase_x_copy:],x_true_copy[:phase_x_copy]))
    y_true_copy = np.concatenate((y_true_copy[phase_y_copy:],y_true_copy[:phase_y_copy]))
    amp_x_copy = (max(x_eqdist_copy)-min(x_eqdist_copy))/2
    amp_y_copy = (max(y_eqdist_copy)-min(y_eqdist_copy))/2
    x_true_copy = amp_x_copy*x_true_copy
    y_true_copy = amp_y_copy*y_true_copy
    # command clocks
    phase_x_command,max_corr_x_command = 0,0
    phase_y_command,max_corr_y_command = 0,0
    for w in range(len(x_true_command)):
        x_shifted = np.concatenate((x_true_command[w:],x_true_command[:w]))
        y_shifted = np.concatenate((y_true_command[w:],y_true_command[:w]))
        corr_x = np.dot(x_eqdist_command,x_shifted)
        corr_y = np.dot(y_eqdist_command,y_shifted)
        if corr_x>max_corr_x_command:
            max_corr_x_command = corr_x
            phase_x_command = w
        if corr_y>max_corr_y_command:
            max_corr_y_command = corr_y
            phase_y_command = w
    #print 'phase for x_true[n] is ',phase_x,' and max correlation is ',max_corr_x
    #print 'phase for y_true[n] is ',phase_y,' and max correlation is ',max_corr_y
    x_true_command = np.concatenate((x_true_command[phase_x_command:],x_true_command[:phase_x_command]))
    y_true_command = np.concatenate((y_true_command[phase_y_command:],y_true_command[:phase_y_command]))
    amp_x_command = (max(x_eqdist_command)-min(x_eqdist_command))/2
    amp_y_command = (max(y_eqdist_command)-min(y_eqdist_command))/2
    x_true_command = amp_x_command*x_true_command
    y_true_command = amp_y_command*y_true_command
    
    # copy clocks
    fig_xt_copy,fig_yt_copy = plt.figure(),plt.figure()
    fig_xt_copy.subplots_adjust(hspace=0.4)
    fig_yt_copy.subplots_adjust(hspace=0.4)
    xt_copy,yt_copy = fig_xt_copy.add_subplot(311),fig_yt_copy.add_subplot(311)
    xt_eqdist_copy,yt_eqdist_copy = fig_xt_copy.add_subplot(312),fig_yt_copy.add_subplot(312)
    xt_copy.plot(x_copy,label='x[n]')
    xt_eqdist_copy.plot(x_eqdist_copy-x_eqdist_copy[0],label='x_eqdist[n]')
    xt_eqdist_copy.plot(x_true_copy-x_true_copy[0],'k-.',lw=3,label='x_true[n]')
    xt_copy.legend(loc='best',frameon=False)
    xt_eqdist_copy.legend(loc='best',frameon=False)
    xt_eqdist_copy.legend(loc='best',frameon=False)
    yt_copy.plot(y_copy,label='y[n]')
    yt_eqdist_copy.plot(y_eqdist_copy-y_eqdist_copy[0],label='y_eqdist[n]')
    yt_eqdist_copy.plot(y_true_copy-y_true_copy[0],'k-.',lw=3,label='y_true[n]')
    yt_copy.legend(loc='best',frameon=False)
    yt_eqdist_copy.legend(loc='best',frameon=False)
    
    dftxk_copy,dftyk_copy = fig_xt_copy.add_subplot(313),fig_yt_copy.add_subplot(313)
    linecolor_dft = 'red'
    dftxk_copy.plot(posfreq_copy,np.abs(dftx_posfreq_copy),linecolor_dft)
    dftyk_copy.plot(posfreq_copy,np.abs(dfty_posfreq_copy),linecolor_dft)
    dftxk_copy.set_ylim(bottom=min(np.abs(dftx_posfreq_copy[1:]))/10,top=max(np.abs(dftx_copy))*10)
    dftyk_copy.set_ylim(bottom=min(np.abs(dfty_posfreq_copy[1:]))/10,top=max(np.abs(dfty_copy))*10)
    dftxk_copy.set_yscale('log')
    dftyk_copy.set_yscale('log')

    # set axis limits
    xt_copy.set_xlim(right=len(x_copy))
    yt_copy.set_xlim(right=len(y_copy))
    xt_eqdist_copy.set_xlim(right=len(x_eqdist_copy))
    yt_eqdist_copy.set_xlim(right=len(y_eqdist_copy))
    dftxk_copy.set_xlim(right=posfreq_copy[-1])
    dftyk_copy.set_xlim(right=posfreq_copy[-1])
    
    # set titles and file labels
    title_fontsize = 20
    #xt_copy.set_title('Sample Copy Clock',fontsize=title_fontsize)
    #yt_copy.set_title('Sample Copy Clock',fontsize=title_fontsize)
    #xy_copy.set_title('Sample Copy Clock',fontsize=title_fontsize)

    fig_xt_copy.text(0.99, 0.96,fname[:len(fname)-4],fontsize=10,color='red',va='baseline',ha='right',multialignment='left')
    fig_yt_copy.text(0.99, 0.96,fname[:len(fname)-4],fontsize=10,color='red',va='baseline',ha='right',multialignment='left')

    if 'YDU' in fname:
        fig_xt_copy.text(0.25, 0.955, 'HEALTHY',fontsize=15,color='black',va='baseline',ha='right',multialignment='left')
        fig_yt_copy.text(0.25, 0.955, 'HEALTHY',fontsize=15,color='black',va='baseline',ha='right',multialignment='left')
    elif 'CIN' in fname:
        fig_xt_copy.text(0.25, 0.955, 'IMPAIRED',fontsize=15,color='black',va='baseline',ha='right',multialignment='left')
        fig_yt_copy.text(0.25, 0.955, 'IMPAIRED',fontsize=15,color='black',va='baseline',ha='right',multialignment='left')
    else:
        print 'not a valid filename'
        
    # set axis labels
    x_axis_fontsize = 15
    xt_copy.set_xlabel('n',fontsize=x_axis_fontsize)
    yt_copy.set_xlabel('n',fontsize=x_axis_fontsize)
    xt_eqdist_copy.set_xlabel('n',fontsize=x_axis_fontsize)
    yt_eqdist_copy.set_xlabel('n',fontsize=x_axis_fontsize)
    dftxk_copy.set_xlabel('k',fontsize=x_axis_fontsize)
    dftyk_copy.set_xlabel('k',fontsize=x_axis_fontsize)

    y_axis_fontsize = 12
    linecolor_xy = 'blue'
    xt_copy.set_ylabel('x[n]',color=linecolor_xy,fontsize=y_axis_fontsize)
    yt_copy.set_ylabel('y[n]',color=linecolor_xy,fontsize=y_axis_fontsize)
    xt_eqdist_copy.set_ylabel('x_eqdist[n],\nx_true[n]',color=linecolor_xy,fontsize=y_axis_fontsize)
    yt_eqdist_copy.set_ylabel('y_eqdist[n],\ny_true[n]',color=linecolor_xy,fontsize=y_axis_fontsize)
    for x1 in xt_copy.get_yticklabels():
        x1.set_color(linecolor_xy)
    for y1 in yt_copy.get_yticklabels():
        y1.set_color(linecolor_xy)
    for x2 in xt_eqdist_copy.get_yticklabels():
        x2.set_color(linecolor_xy)
    for y2 in yt_eqdist_copy.get_yticklabels():
        y2.set_color(linecolor_xy)
        
    dftxk_copy.set_ylabel(r'$|X[k]|$',color=linecolor_dft,fontsize=y_axis_fontsize)
    dftyk_copy.set_ylabel(r'$|Y[k]|$',color=linecolor_dft,fontsize=y_axis_fontsize)
    for v1 in dftxk_copy.get_yticklabels():
        v1.set_color(linecolor_dft)
    for v2 in dftyk_copy.get_yticklabels():
        v2.set_color(linecolor_dft)
    
    # save figures
    fig_xt_copy.savefig(path+'figs_raw/'+fname[:len(fname)-4]+'/x_true_copy_'+fname[:len(fname)-4]+'.png')
    fig_yt_copy.savefig(path+'figs_raw/'+fname[:len(fname)-4]+'/y_true_copy_'+fname[:len(fname)-4]+'.png')
    
    # command clocks
    fig_xt_command,fig_yt_command = plt.figure(),plt.figure()
    fig_xt_command.subplots_adjust(hspace=0.4)
    fig_yt_command.subplots_adjust(hspace=0.4)
    xt_command,yt_command = fig_xt_command.add_subplot(311),fig_yt_command.add_subplot(311)
    xt_eqdist_command,yt_eqdist_command = fig_xt_command.add_subplot(312),fig_yt_command.add_subplot(312)
    xt_command.plot(x_command,label='x[n]')
    xt_eqdist_command.plot(x_eqdist_command-x_eqdist_command[0],label='x_eqdist[n]')
    xt_eqdist_command.plot(x_true_command-x_true_command[0],'k-.',lw=3,label='x_true[n]')
    xt_command.legend(loc='best',frameon=False)
    xt_eqdist_command.legend(loc='best',frameon=False)
    xt_eqdist_command.legend(loc='best',frameon=False)
    yt_command.plot(y_command,label='y[n]')
    yt_eqdist_command.plot(y_eqdist_command-y_eqdist_command[0],label='y_eqdist[n]')
    yt_eqdist_command.plot(y_true_command-y_true_command[0],'k-.',lw=3,label='y_true[n]')
    yt_command.legend(loc='best',frameon=False)
    yt_eqdist_command.legend(loc='best',frameon=False)
    
    dftxk_command,dftyk_command = fig_xt_command.add_subplot(313),fig_yt_command.add_subplot(313)
    linecolor_dft = 'red'
    dftxk_command.plot(posfreq_command,np.abs(dftx_posfreq_command),linecolor_dft)
    dftyk_command.plot(posfreq_command,np.abs(dfty_posfreq_command),linecolor_dft)
    dftxk_command.set_ylim(bottom=min(np.abs(dftx_posfreq_command[1:]))/10,top=max(np.abs(dftx_command))*10)
    dftyk_command.set_ylim(bottom=min(np.abs(dfty_posfreq_command[1:]))/10,top=max(np.abs(dfty_command))*10)
    dftxk_command.set_yscale('log')
    dftyk_command.set_yscale('log')

    # set axis limits
    xt_command.set_xlim(right=len(x_command))
    yt_command.set_xlim(right=len(y_command))
    xt_eqdist_command.set_xlim(right=len(x_eqdist_command))
    yt_eqdist_command.set_xlim(right=len(y_eqdist_command))
    dftxk_command.set_xlim(right=posfreq_command[-1])
    dftyk_command.set_xlim(right=posfreq_command[-1])
    
    # set titles and file labels
    title_fontsize = 20
    #xt_command.set_title('Sample command Clock',fontsize=title_fontsize)
    #yt_command.set_title('Sample command Clock',fontsize=title_fontsize)
    #xy_command.set_title('Sample command Clock',fontsize=title_fontsize)

    fig_xt_command.text(0.99, 0.96,fname[:len(fname)-4],fontsize=10,color='red',va='baseline',ha='right',multialignment='left')
    fig_yt_command.text(0.99, 0.96,fname[:len(fname)-4],fontsize=10,color='red',va='baseline',ha='right',multialignment='left')

    if 'YDU' in fname:
        fig_xt_command.text(0.25, 0.955, 'HEALTHY',fontsize=15,color='black',va='baseline',ha='right',multialignment='left')
        fig_yt_command.text(0.25, 0.955, 'HEALTHY',fontsize=15,color='black',va='baseline',ha='right',multialignment='left')
    elif 'CIN' in fname:
        fig_xt_command.text(0.25, 0.955, 'IMPAIRED',fontsize=15,color='black',va='baseline',ha='right',multialignment='left')
        fig_yt_command.text(0.25, 0.955, 'IMPAIRED',fontsize=15,color='black',va='baseline',ha='right',multialignment='left')
    else:
        print 'not a valid filename'
        
    # set axis labels
    x_axis_fontsize = 15
    xt_command.set_xlabel('n',fontsize=x_axis_fontsize)
    yt_command.set_xlabel('n',fontsize=x_axis_fontsize)
    xt_eqdist_command.set_xlabel('n',fontsize=x_axis_fontsize)
    yt_eqdist_command.set_xlabel('n',fontsize=x_axis_fontsize)
    dftxk_command.set_xlabel('k',fontsize=x_axis_fontsize)
    dftyk_command.set_xlabel('k',fontsize=x_axis_fontsize)

    y_axis_fontsize = 12
    linecolor_xy = 'blue'
    xt_command.set_ylabel('x[n]',color=linecolor_xy,fontsize=y_axis_fontsize)
    yt_command.set_ylabel('y[n]',color=linecolor_xy,fontsize=y_axis_fontsize)
    xt_eqdist_command.set_ylabel('x_eqdist[n],\nx_true[n]',color=linecolor_xy,fontsize=y_axis_fontsize)
    yt_eqdist_command.set_ylabel('y_eqdist[n],\ny_true[n]',color=linecolor_xy,fontsize=y_axis_fontsize)
    for x1 in xt_command.get_yticklabels():
        x1.set_color(linecolor_xy)
    for y1 in yt_command.get_yticklabels():
        y1.set_color(linecolor_xy)
    for x2 in xt_eqdist_command.get_yticklabels():
        x2.set_color(linecolor_xy)
    for y2 in yt_eqdist_command.get_yticklabels():
        y2.set_color(linecolor_xy)
        
    dftxk_command.set_ylabel(r'$|X[k]|$',color=linecolor_dft,fontsize=y_axis_fontsize)
    dftyk_command.set_ylabel(r'$|Y[k]|$',color=linecolor_dft,fontsize=y_axis_fontsize)
    for v1 in dftxk_command.get_yticklabels():
        v1.set_color(linecolor_dft)
    for v2 in dftyk_command.get_yticklabels():
        v2.set_color(linecolor_dft)

    # save figures
    fig_xt_command.savefig(path+'figs_raw/'+fname[:len(fname)-4]+'/x_true_command_'+fname[:len(fname)-4]+'.png')
    fig_yt_command.savefig(path+'figs_raw/'+fname[:len(fname)-4]+'/y_true_command_'+fname[:len(fname)-4]+'.png')

    plt.close('all')

# compare percent energy in peak for the drawings of healthy vs. impaired patients
#binedges = ct.get_bins(energy_in_peak_copy,nbins=10)
#ct.make_hist([elt[1] for elt in energy_in_peak_copy if elt[0]=='healthy'],
             #[elt[1] for elt in energy_in_peak_copy if elt[0]=='impaired'],
             #binedges,'Percent Energy in Largest DFS Coefficient','energy_in_peak_copy',path)
#binedges = ct.get_bins(energy_in_peak_command,nbins=10)
#ct.make_hist([elt[1] for elt in energy_in_peak_command if elt[0]=='healthy'],
             #[elt[1] for elt in energy_in_peak_command if elt[0]=='impaired'],
             #binedges,'Percent Energy in Largest DFS Coefficient','energy_in_peak_command',path)
