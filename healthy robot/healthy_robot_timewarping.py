# for sanity checks on normalizing for non-constant velocity

# focus on comparing one healthy robot's drawing to the
# ideal underlying sinusoid (so should choose N = size of
# DFT = N_orig, where N_orig is the number of points in
# original drawing)

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

def sinc(omega_c,n,length_of_sinc):
    if n==0:
        return 1
    elif abs(n)>(length_of_sinc-1)/2:
        return 0
    else:
        return math.sin(omega_c*n)/(omega_c*n)

# get coordinates and timestamps
n = np.arange(142)
noise_x = np.random.normal(loc=0,scale=0.01,size=len(n))
noise_y = np.random.normal(loc=0,scale=0.01,size=len(n))
x = np.cos(1./2.*(2.*pi*n/250.)**2)+noise_x # for a healthy robot, comment out the noise
y = np.sin(1./2.*(2.*pi*n/250.)**2)+noise_y
# make the starting point was something other than (1,0)
offset = 100
x = np.concatenate((x[offset:],x[:offset]))
y = np.concatenate((y[offset:],y[:offset]))
N_orig = len(x)

# velocity/distance (T = 1)
dists = []
for w in range(1,len(x)):
    dx,dy = x[w]-x[w-1],y[w]-y[w-1]
    dist = math.sqrt(dx**2+dy**2)
    dists.append(dist)
dist_avg = mean(dists)
print 'average distance between points is ',dist_avg
print 'total distance is ',sum(dists)

# want to get 142 evenly-spaced points along the curve

# generate a much longer array with 199 linearly-interpolated 
# points between the actual data points
x_interp,y_interp = [],[]
for w in range(len(x)-1):
    x_interp.append(x[w])
    y_interp.append(y[w])
    dx,dy = x[w+1]-x[w],y[w+1]-y[w]
    dist = math.sqrt(dx**2+dy**2)
    n_segments = ceil(dist/dist_avg)*100
    for r in range(1,int(n_segments)):
        x_new = x[w]+r*dx/n_segments
        y_new = y[w]+r*dy/n_segments
        x_interp.append(x_new)
        y_interp.append(y_new)
x_interp.append(x[-1])
y_interp.append(y[-1])
# check
dists_interp = []
for w in range(1,len(x_interp)):
    dx,dy = x_interp[w]-x_interp[w-1],y_interp[w]-y_interp[w-1]
    dist = math.sqrt(dx**2+dy**2)
    dists_interp.append(dist)
dist_avg_interp = mean(dists_interp)
print '\naverage distance between interpolated points is ',dist_avg_interp
print 'total distance is now ',sum(dists_interp)

# start from the first point and find the ones that are 
# approximately a distance dist_avg from each other
x_eqdist,y_eqdist = [x_interp[0]],[y_interp[0]]
idx = 0
for k in range(len(x)-1):
    dist_total = 0
    for j in range(idx,len(x_interp)-1):
        dx,dy = x_interp[j+1]-x_interp[j],y_interp[j+1]-y_interp[j]
        dist_total += math.sqrt(dx**2+dy**2)
        #if abs(dist_total-dist_avg)<0.0005:
        if abs(dist_total-dist_avg)<dist_avg/100.:
            idx = j+1
            break
    x_eqdist.append(x_interp[idx])
    y_eqdist.append(y_interp[idx])
x_eqdist,y_eqdist = np.array(x_eqdist),np.array(y_eqdist)
# check
dists_check = []
for w in range(1,len(x_eqdist)):
    dx,dy = x_eqdist[w]-x_eqdist[w-1],y_eqdist[w]-y_eqdist[w-1]
    dist = math.sqrt(dx**2+dy**2)
    dists_check.append(dist)
dist_avg_check = mean(dists_check)
print '\naverage distance between points is now ',dist_avg_check
print 'total distance is now ',sum(dists_check)

# now want to estimate the frequency of the underlying sinusoid
# subtract mean values (zm = zero mean)
x_eqdist_zm,y_eqdist_zm = x_eqdist-mean(x_eqdist),y_eqdist-mean(y_eqdist)
    
# DFS coefficients
dft_size = N_orig
k = range(dft_size)
dftx,dfty = np.fft.fft(x_eqdist_zm,n=dft_size),np.fft.fft(y_eqdist_zm,n=dft_size)

# k_near_pi is the smallest k value for which w_k = 2*pi*k/N is
# greater than pi
k_near_pi = 0
if dft_size%2==0:
    k_near_pi = dft_size/2+1
else:
    k_near_pi = math.ceil(dft_size/2)

# only use the positive frequencies for plotting
pos_k = k[:int(k_near_pi)]
dftx_pos_k,dfty_pos_k = dftx[:k_near_pi],dfty[:k_near_pi]

# take the frequency of the largest DFS coefficient to be
# the approximate frequency of the underlying sinusoid
abs_dftx_pos_k,abs_dfty_pos_k = list(np.abs(dftx_pos_k)),list(np.abs(dfty_pos_k))
k_true_x,k_true_y = abs_dftx_pos_k.index(max(abs_dftx_pos_k)),abs_dfty_pos_k.index(max(abs_dfty_pos_k))
w_true_x,w_true_y = 2*pi*k_true_x/dft_size,2*pi*k_true_y/dft_size
print '\ncheck: ',w_true_x,' = ',2.*pi/dft_size,' = ',w_true_y,'?'
print '\n',np.angle(dftx[1]),np.angle(dftx[-1]),' ---> ',np.angle(dftx[1])-np.angle(dftx[-1])
print np.angle(dfty[1]),np.angle(dfty[-1]),' ---> ',np.angle(dfty[1])-np.angle(dfty[-1])

# calculate ideal underlying sinusoid
x_true = np.cos(w_true_x*n)
y_true = np.sin(w_true_y*n)
# use maximum correlation to determine phase
phase_x,max_corr_x = 0,0
phase_y,max_corr_y = 0,0
for w in range(len(x_true)):
    x_shifted = np.concatenate((x_true[w:],x_true[:w]))
    y_shifted = np.concatenate((y_true[w:],y_true[:w]))
    corr_x = np.dot(x_eqdist,x_shifted)
    corr_y = np.dot(y_eqdist,y_shifted)
    if corr_x>max_corr_x:
        max_corr_x = corr_x
        phase_x = w
    if corr_y>max_corr_y:
        max_corr_y = corr_y
        phase_y = w
print 'phase for x_true[n] is ',phase_x,' and max correlation is ',max_corr_x
print 'phase for y_true[n] is ',phase_y,' and max correlation is ',max_corr_y
x_true = np.concatenate((x_true[phase_x:],x_true[:phase_x]))
y_true = np.concatenate((y_true[phase_y:],y_true[:phase_y]))

# compute correlation between the two
corr_x,corr_y = np.dot(x_eqdist,x_true),np.dot(y_eqdist,y_true)
Etrue_x,Etrue_y = sum(x_true**2),sum(y_true**2)
print '\n correlation for x[n] = ',corr_x,' and correlation for y[n] = ',corr_y
print '\n energy in x_true[n] = ',Etrue_x,' and energy in y_true[n] = ',Etrue_y
# compute mean squared difference between the two
Ediff_x,Ediff_y = sum((x_eqdist-x_true)**2),sum((y_eqdist-y_true)**2)
print '\n mean square difference for x[n] = ',Ediff_x,' and mean square difference for y[n] = ',Ediff_y

plt.close('all')

fig1 = plt.figure()
ax1 = fig1.add_subplot(211)
ax1.plot(x,label='x')
ax1.plot(y,label='y')
ax1.legend(loc='best',frameon=False)
ax1.set_ylabel('x,y',fontsize=20)
ax1.set_xlim(right=len(x))
ax1.set_ylim(bottom=min(x)-0.2,top=max(x)+0.2)
ax2 = fig1.add_subplot(212)
ax2.plot(dists)
ax2.set_ylabel('v',fontsize=20)
fig1.savefig(path+'figs_raw/healthy_robot/xy_noisy_nonconstant_velocity.png')

fig2 = plt.figure()
ax3 = fig2.add_subplot(211)
ax3.plot(x_eqdist,label='x_eqdist')
ax3.plot(y_eqdist,label='y_eqdist')
ax3.legend(loc='best',frameon=False)
ax3.set_ylabel('x_eqdist,\ny_eqdist',fontsize=15)
ax3.set_xlim(right=len(x_eqdist))
ax4 = fig2.add_subplot(212)
ax4.plot(dists_check)
ax4.set_ylabel('v_eqdist',fontsize=20)
ax4.set_ylim(bottom=min(dists),top=max(dists))
fig2.savefig(path+'figs_raw/healthy_robot/xy_noisy_constant_velocity.png')

fig3 = plt.figure()
ax5 = fig3.add_subplot(211)
ax5.stem(pos_k,np.abs(dftx_pos_k))
ax5.set_xlabel('k',fontsize=20)
ax5.set_ylabel('|X[k]|',fontsize=20)
ax5.set_xlim(right=max(pos_k)+1)
ax6 = fig3.add_subplot(212)
ax6.stem(pos_k,np.abs(dfty_pos_k))
ax6.set_xlabel('k',fontsize=20)
ax6.set_ylabel('|Y[k]|',fontsize=20)
ax6.set_xlim(right=max(pos_k)+1)
fig3.savefig(path+'figs_raw/healthy_robot/xy_noisy_constant_velocity_DFS.png')

fig4 = plt.figure()
ax7 = fig4.add_subplot(211)
ax7.plot(x_eqdist,label='x_eqdist')
ax7.plot(x_true,label='x_true')
ax7.legend(loc='best',frameon=False)
ax7.set_ylabel('x_eqdist,\nx_true',fontsize=15)
ax7.set_xlim(right=len(x))
ax8 = fig4.add_subplot(212)
ax8.plot(y_eqdist,label='y_eqdist')
ax8.plot(y_true,label='y_true')
ax8.legend(loc='best',frameon=False)
ax8.set_xlabel('n',fontsize=20)
ax8.set_ylabel('y_eqdist,\ny_true',fontsize=15)
ax8.set_xlim(right=len(y))
fig4.savefig(path+'figs_raw/healthy_robot/xy_noisy_true_sinusoid.png')

fig_xy = plt.figure()
ax_xy = fig_xy.add_subplot(111)
ax_xy.plot(x_eqdist,y_eqdist,label='x_eqdist,y_eqdist')
ax_xy.legend(loc='best',frameon=False)
ax_xy.set_xlabel('x_eqdist')
ax_xy.set_ylabel('y_eqdist')
fig_xy.savefig(path+'figs_raw/healthy_robot/circle_noisy_nonconstant_velocity.png')

plt.show()
