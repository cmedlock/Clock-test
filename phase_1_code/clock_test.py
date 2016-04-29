import math
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as anim
import numpy as np
import numpy.fft
import os
from pylab import *

def parse_file(fname,path):
    # input: data file
    # output: x,y,p,t,
    # x_command,y_command,p_command,t_command (lists)

    print 'reading file ',fname,'...'

    f = open(path+fname)
    data = f.readlines()

    # copy or command clock?
    clock_type = ''
    x,y,p,t = [],[],[],[]
    x_command,y_command,p_command,t_command = [],[],[],[]

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
            if 'symbol label' in line and clock_type=='COPY' and len(x)>0:
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
        pressure = double(line[5])
        timestamp = double(line[7])
        if clock_type=='COPY':
            # if a point has been skipped, leave a space for it with the appropriate timestamp
            if len(x)>0 and timestamp-t[-1]>20:
                x.append(-5)
                y.append(-5)
                p.append(-15)
                t.append(t[-1]+(timestamp-t[-1])/2)
            x.append(xcoord)
            y.append(ycoord)
            p.append(pressure)
            t.append(timestamp)
        elif clock_type=='COMMAND':
            if len(x_command)>0 and timestamp-t_command[-1]>20:
                x_command.append(-5)
                y_command.append(-5)
                p_command.append(-5)
                t_command.append(t_command[-1]+(timestamp-t_command[-1])/2)
            x_command.append(xcoord)
            y_command.append(ycoord)
            p_command.append(pressure)
            t_command.append(timestamp)
        else:
            print 'not a valid clock type'
    
    f.close()

    return x,y,p,t,x_command,y_command,p_command,t_command

def interpolate(x):
    # input: x (list) with missing points marked by the value -5
    # output: xintp (lists) with the missing points filled in
    #         using a sinc function for interpolation

    if -5 in x:
        missing_points = []
        for w in range(len(x)):
            if x[w]==-5:
                missing_points.append(w)
        #print '   list x is missing points ',missing_points,', downsampling and expanding...'
        # interpolate to find the missing points
        xintp = [elt for elt in x]
        for missing_point_init in missing_points:
            #print '   looking for point ',missing_point_init
            sinc = []
            # expand by 2 (insert zeros)
            # downsample then expand by 2
            # this code can only handle a single missing point with an odd index
            expanded = []
            # if the missing point had an even index, want to save all of the odd samples
            # do this by appending a zero to the start of the list before downsampling and
            # expanding it
            missing_point = missing_point_init
            if missing_point_init%2==0:
                x_shiftedbyone = [0]+list(x)
                missing_point_oddidx = x_shiftedbyone.index(-5)
                #print '   missing point has been shifted to index ',missing_point_oddidx
                missing_point = missing_point_oddidx
                for d in range(len(x_shiftedbyone)):
                    xval = x_shiftedbyone[d]
                    if d%2==0:
                        expanded.append(xval if xval!=-5 else 0)
                    else:
                        expanded.append(0)
            # if the missing point had an odd index, want to save all of the even samples
            elif missing_point_init%2==1:
                for d in range(len(x)):
                    xval = x[d]
                    if d%2==0:
                        expanded.append(xval if xval!=-5 else 0)
                    else:
                        expanded.append(0)
            #print x[1:10]
            #print expanded[:10]
            # windowed sinc filter
            L = 2
            n_neg = np.arange(-missing_point,0)
            n_pos = np.arange(1,len(expanded)-missing_point)
            sinc_neg = np.sin(np.pi*n_neg/L)/(np.pi*n_neg/L)
            sinc_pos = np.sin(np.pi*n_pos/L)/(np.pi*n_pos/L)
            # NB: normally would set sinc[0] equal to 1, but in this
            # case we don't want the -5 at that point to contribute to
            # the sum, so just set sinc[0] equal to 0
            sinc = np.concatenate((sinc_neg,np.array([0]),sinc_pos))
            # evaluate convolution at missing point
            #print '*** ',expanded.count(-5)
            missing = np.dot(expanded,sinc)
            xintp[missing_point_init] = missing
        return xintp
    return x

def linear_interpolate(x):
    # input: x (list) with missing points marked by the value -5
    # output: xintp (list) with the missing points filled in
    #         using linear interpolation

    if -5 in x:
        xintp = []
        for w in range(len(x)):
            if x[w]!=-5:
                xintp.append(x[w])
            else:
                xintp.append((x[w-1]+x[w+1])/2)
        return xintp
    return x
        
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

def draw_clock(x,y,xname,yname,name,fname,path):
    # input: x[t],y[t],t,x- and y-axis labels
    #        fname,path in order to save the figure and video
    # output: xy_name_fname.png with plot of circle

    # figures
    fig_xy = plt.figure()
    xy = fig_xy.add_subplot(111)
    xy.plot(x,y)

    # equalize axis scales for the clock drawings (not the pressure plots)
    if xname=='x' and yname=='y':
        xy.set_aspect(1.0)
        if max(x)-min(x)>max(y)-min(y):
            ax_range = max(x)-min(x)+20
            xy.set_xlim(min(x)-10,max(x)+10)
            xy.set_ylim(min(y)-10,min(y)+ax_range)
        else:
            ax_range = max(y)-min(y)+20
            xy.set_xlim(min(x)-10,min(x)+ax_range)
            xy.set_ylim(min(y)-10,max(y)+10)
        plt.axis('equal')
    
    # set titles and file labels
    title_fontsize = 20
    #xt.set_title('Sample Copy Clock',fontsize=title_fontsize)
    #yt.set_title('Sample Copy Clock',fontsize=title_fontsize)
    #xy.set_title('Sample Copy Clock',fontsize=title_fontsize)

    fig_xy.text(0.99, 0.96, fname[:len(fname)-4],fontsize=10,color='red',va='baseline',ha='right',multialignment='left')

    if 'YDU' in fname:
        fig_xy.text(0.25, 0.955, 'HEALTHY',fontsize=15,color='black',va='baseline',ha='right',multialignment='left')
    elif 'CIN' in fname:
        fig_xy.text(0.25, 0.955, 'IMPAIRED',fontsize=15,color='black',va='baseline',ha='right',multialignment='left')
    else:
        print 'not a valid filename'

    # set axis labels
    x_axis_fontsize = 20
    xy.set_xlabel(xname,fontsize=x_axis_fontsize)

    y_axis_fontsize = 20
    xy.set_ylabel(yname,fontsize=y_axis_fontsize)

    # save figures
    fig_xy.savefig(path+'figs_raw/'+fname[:len(fname)-4]+'/'+name+'_'+fname[:len(fname)-4]+'.png')
                      
def plot_xyt_other(x,t,xname,tname,otherx,w,otherxname,wname,logplot,othername,fname,path):
    # input: x[t],t,x- and y-axis labels
    #        other[w],w (for e.g. |X[k]|,k), x- and y-axis labels
    #        log plot option for spectra
    #        othername,fname,path in order to save the figure
    # output: xt_othername_fname.png with plot of x[t]
    #         directly above plot of other[n]

    fig_xt = plt.figure()
    fig_xt.subplots_adjust(left=0.15,hspace=0.3)
    xt = fig_xt.add_subplot(211)
    xt.plot(t,x)
    
    otherxw = fig_xt.add_subplot(212)  
    linecolor_other = 'red'
    otherxw.plot(w,otherx,linecolor_other)
    if logplot:
        otherxw.set_ylim(bottom=min(otherx[1:])/10,top=max(otherx)*10)
        otherxw.set_yscale('log')

    # set titles and file labels
    title_fontsize = 20
    #xt.set_title('Sample Copy Clock',fontsize=title_fontsize)

    fig_xt.text(0.99, 0.96,fname[:len(fname)-4],fontsize=10,color='red',va='baseline',ha='right',multialignment='left')

    if 'YDU' in fname:
        fig_xt.text(0.25, 0.955, 'HEALTHY',fontsize=15,color='black',va='baseline',ha='right',multialignment='left')
    elif 'CIN' in fname:
        fig_xt.text(0.25, 0.955, 'IMPAIRED',fontsize=15,color='black',va='baseline',ha='right',multialignment='left')
    else:
        print 'not a valid filename'
        
    # set axis labels
    x_axis_fontsize = 20
    xt.set_xlabel('time [ms]',fontsize=x_axis_fontsize)
    otherxw.set_xlabel(wname,fontsize=x_axis_fontsize)

    y_axis_fontsize = 25
    linecolor_xy = 'blue'
    xt.set_ylabel(xname,color=linecolor_xy,fontsize=y_axis_fontsize)
    for x1 in xt.get_yticklabels():
        x1.set_color(linecolor_xy)
        
    otherxw.set_ylabel(r'$'+otherxname+'$',color=linecolor_other,fontsize=y_axis_fontsize)
    for d1 in otherxw.get_yticklabels():
        d1.set_color(linecolor_other)
    
    # save figures
    #fig_xt.savefig(path+'figs_raw/'+fname[:len(fname)-4]+'/'+othername+'_'+fname[:len(fname)-4]+'.png')
    fig_xt.savefig(path+'figs_raw/'+fname+'/'+othername+'_'+fname+'.png')
    
    plt.close('all')

def find_relsidelobe_amplitude(dft):
    # input: dft (list) = |X[k]| for the positive frequences for a real signal x[n]
    # output: rel_sidelobe_amplitude = relative side lobe amplitude of dft in dB
    abs_dft = np.abs(dft)[:len(dft)/2] # restrict to only the positive frequencies
    peak_amplitude,sidelobe_amplitude = 0,0
    for w in range(1,len(abs_dft)-1):
        if abs_dft[w]-abs_dft[w-1]>0 and abs_dft[w]-abs_dft[w+1]>0:
            if peak_amplitude==0:
                peak_amplitude = abs_dft[w]
            elif peak_amplitude!=0 and sidelobe_amplitude==0:
                sidelobe_amplitude = abs_dft[w]
        if peak_amplitude!=0 and sidelobe_amplitude!=0:
            break
    relsidelobe_amplitude = 20*math.log10(peak_amplitude/sidelobe_amplitude)
    return relsidelobe_amplitude

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
    #n,bins,patches = hist.hist([healthy,impaired],nbins,histtype='bar',label=['Healthy','Impaired'])
    #n,bins_,patches = hist.hist([healthy,impaired],bins=10 ** np.linspace(np.log10(range_min), np.log10(range_max), 10),histtype='bar',label=['Healthy','Impaired'])
    n,bins,patches = hist.hist([healthy,impaired],nbins,range=[range_min,range_max],histtype='bar',label=['Healthy','Impaired'])
    log_scale_hist = np.log10(n)
    print len(log_scale_hist)
    fig_2 = plt.figure()
    hist_2 = fig_2.add_subplot(111)
    n_,bins,patches = hist_2.hist(log_scale_hist,nbins,range=[range_min,range_max],histtype='bar',label=['Healthy','Impaired'])
    n1,n2 = n[0],n[1]
    hist_2.set_ylim(top=max(max(n1),max(n2))*1.2)
    hist_2.legend(loc='best',frameon=False)
    
    # set titles and file labels
    title_fontsize = 15
    #hist.set_title(nameforplot+': Healthy vs. Impaired',fontsize=title_fontsize)
    #hist.set_title(nameforplot,fontsize=title_fontsize)

    # set axis label
    x_axis_fontsize = 15
    hist_2.set_xlabel(nameforplot,fontsize=x_axis_fontsize)

    # save figure
    fig_2.savefig(path+'/compare_healthy_impaired/compare_'+nameforfile+'.png')
    
    plt.close('all')
    