import math
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import numpy.fft
import os
from pylab import *

def parse_file(path,fname):
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

def interpolate(x,y):
    # input: x,y (lists) with missing points marked by the value -5
    # output: xintp,yintp (lists) with the missing points filled in

    if len(x)!=len(y):
        raise NameError('ERROR: x and y lists do not have the same length')
    
    sinc = None
    xexpand,yexpand = [],[]
    if -5 in x:
        # expand by 2 (insert zeros)
        print '   lists x and y have missing points, expanding a[] and b[t]...'
        for v in range(len(x)-1):
            xexpand.append(x[v])
            yexpand.append(y[v])
            if x[v]!=-5 and x[v+1]!=-5:
                xexpand.append(0)
                yexpand.append(0)
        xexpand.append(x[-1])
        yexpand.append(y[-1])
        for d in range(len(xexpand)):
            if xexpand[d]!=-5:
                continue
            #print '      found missing point at index ',(d-1)/2
            # windowed sinc filter
            #print '         forming windowed sinc filter...'
            L = 2
            n_neg = np.arange(-d,0)
            n_pos = np.arange(1,len(xexpand)-d)
            sinc_neg = np.sin(np.pi*n_neg/L)/(np.pi*n_neg/L)
            sinc_pos = np.sin(np.pi*n_pos/L)/(np.pi*n_pos/L)
            # NB: normally would set sinc[0] equal to 1, but in this
            # case we don't want the -5 at that point to contribute to
            # the sum, so just set sinc[0] equal to 0
            sinc = np.concatenate((sinc_neg,np.array([0]),sinc_pos))
            # evaluate convolution at missing point
            missing_x = np.dot(xexpand,sinc)
            missing_y = np.dot(yexpand,sinc)
            #print '         missing x is ',missing_x,' and missing y is ',missing_y
            xexpand[d] = missing_x
            yexpand[d] = missing_y
        # downsample by 2 to remove the zeros
        xintp,yintp = [],[]
        for v in range(len(xexpand)):
            if v%2==0:
                xintp.append(xexpand[v])
                yintp.append(yexpand[v])
            elif v%2==1 and xexpand[v]!=0:
                #print '      inserting missing point (',xexpand[v],',',yexpand[v],') at index ',v/2
                xintp.append(xexpand[v])
                yintp.append(yexpand[v])
        return xintp,yintp
    return x,y

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
        ax.append(a_y)
    return vx,vy,ax,ay
      
def plot_xyt_other(x,y,t,otherx,othery,n,othername,fname):
    # input: x[t],y[t],t
    #        other[n],n (for e.g. |X[k]|,k)
    # output: xyt_othername_fname.png with plots of x[t] or y[t]
    #         directly above other[n]

    fig_xt,fig_yt = plt.figure(),plt.figure()
    fig_xt.subplots_adjust(right=0.85)
    fig_yt.subplots_adjust(right=0.85)

    xt,yt = fig_xt.add_subplot(211),fig_yt.add_subplot(211)
    xt.plot(t-t[0],x)
    yt.plot(t-t[0],y)

    vxt,vyt = fig_xt.add_subplot(212),fig_yt.add_subplot(212)    
    marker_size,marker_edgecolor_v,linecolor_v = 5,'darkred','red'
    vxt.plot(t[1:]-t[0],vx,linecolor_v,markersize=marker_size)#,markeredgecolor=marker_edgecolor_v)
    vyt.plot(t[1:]-t[0],vy,linecolor_v,markersize=marker_size)#,markeredgecolor=marker_edgecolor_v)

    axt,ayt = vxt.twinx(),vyt.twinx()
    marker_edgecolor_a,linecolor_a = 'darkgreen','green'
    axt.plot(t[2:]-t[0],ax,linecolor_a,markersize=marker_size)#,markeredgecolor=marker_edgecolor_a)
    ayt.plot(t[2:]-t[0],ay,linecolor_a,markersize=marker_size)#,markeredgecolor=marker_edgecolor_a)
    
    # set titles and file labels
    title_fontsize = 20
    #xt.set_title('Sample Copy Clock',fontsize=title_fontsize)
    #yt.set_title('Sample Copy Clock',fontsize=title_fontsize)
    #xy.set_title('Sample Copy Clock',fontsize=title_fontsize)

    fig_xt.text(0.99, 0.96,fname[:len(fname)-4],fontsize=10,color='red',va='baseline',ha='right',multialignment='left')
    fig_yt.text(0.99, 0.96,fname[:len(fname)-4],fontsize=10,color='red',va='baseline',ha='right',multialignment='left')

    if 'YDU' in fname:
        fig_xt.text(0.25, 0.955, 'HEALTHY',fontsize=15,color='black',va='baseline',ha='right',multialignment='left')
        fig_yt.text(0.25, 0.955, 'HEALTHY',fontsize=15,color='black',va='baseline',ha='right',multialignment='left')
    elif 'CIN' in fname:
        fig_xt.text(0.25, 0.955, 'IMPAIRED',fontsize=15,color='black',va='baseline',ha='right',multialignment='left')
        fig_yt.text(0.25, 0.955, 'IMPAIRED',fontsize=15,color='black',va='baseline',ha='right',multialignment='left')
    else:
        print 'not a valid filename'
        
    # set axis labels
    x_axis_fontsize = 20
    vxt.set_xlabel('time [ms]',fontsize=x_axis_fontsize)
    vyt.set_xlabel('time [ms]',fontsize=x_axis_fontsize)

    y_axis_fontsize = 25
    linecolor_xy = 'blue'
    xt.set_ylabel(r'$x$',color=linecolor_xy,fontsize=y_axis_fontsize)
    yt.set_ylabel(r'$y$',color=linecolor_xy,fontsize=y_axis_fontsize)
    for x1 in xt.get_yticklabels():
        x1.set_color(linecolor_xy)
    for y1 in yt.get_yticklabels():
        y1.set_color(linecolor_xy)
        
    vxt.set_ylabel(r'$v_x$',color=linecolor_v,fontsize=y_axis_fontsize)
    vyt.set_ylabel(r'$v_y$',color=linecolor_v,fontsize=y_axis_fontsize)
    axt.set_ylabel(r'$a_x$',color=linecolor_a,fontsize=y_axis_fontsize)
    ayt.set_ylabel(r'$a_y$',color=linecolor_a,fontsize=y_axis_fontsize)
    for v1 in vxt.get_yticklabels():
        v1.set_color(linecolor_v)
    for v2 in vyt.get_yticklabels():
        v2.set_color(linecolor_v)
    for m1 in axt.get_yticklabels():
        m1.set_color(linecolor_a)
    for m2 in ayt.get_yticklabels():
        m2.set_color(linecolor_a)
    
    # save figures
    fig_xt.savefig(path+'figs_raw/'+fname[:len(fname)-4]+'/vxt_axt_'+fname[:len(fname)-4]+'.png')
    fig_yt.savefig(path+'figs_raw/'+fname[:len(fname)-4]+'/vyt_ayt_'+fname[:len(fname)-4]+'.png')
    