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

def interpolate(x,y):
    # input: x,y (lists) with missing points marked by the value -5
    # output: xintp,yintp (lists) with the missing points filled in
    #         using a sinc function for interpolation

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
                print '      inserting missing point (',xexpand[v],',',yexpand[v],') at index ',v/2
                xintp.append(xexpand[v])
                yintp.append(yexpand[v])
        return xintp,yintp
    return x,y

def linear_interpolate(x,y):
    # input: x,y (lists) with missing points marked by the value -5
    # output: xintp,yintp (lists) with the missing points filled in
    #         using linear interpolation

    if len(x)!=len(y):
        raise NameError('ERROR: x and y lists do not have the same length')

    if -5 in x:
        xintp,yintp = [],[]
        for w in range(len(x)):
            if x[w]!=-5:
                xintp.append(x[w])
                yintp.append(y[w])
            else:
                xintp.append((x[w-1]+x[w+1])/2)
                yintp.append((y[w-1]+y[w+1])/2)
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
        ay.append(a_y)
    return vx,vy,ax,ay

def draw_clock(x,y,name,fname,path):
    # input: x[t],y[t],t,x- and y-axis labels
    #        fname,path in order to save the figure and video
    # output: xy_name_fname.png with plot of circle,

    # figures
    fig_xy = plt.figure()
    xy = fig_xy.add_subplot(111,aspect=1.0)
    xy.plot(x,y)

    # equalize axis scales
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
    xy.set_xlabel('x',fontsize=x_axis_fontsize)

    y_axis_fontsize = 20
    xy.set_ylabel('y',fontsize=y_axis_fontsize)

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
        otherxw.set_ylim(bottom=min(other[1:])/10,top=max(other)*10)
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
    fig_xt.savefig(path+'figs_raw/'+fname[:len(fname)-4]+'/'+othername+'_'+fname[:len(fname)-4]+'.png')
    
    plt.close('all')
    