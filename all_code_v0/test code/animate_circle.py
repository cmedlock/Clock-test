import math
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as anim
import numpy as np
import os
from pylab import *

plt.close('all')

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xlim(-2,2)
ax.set_ylim(-2,2)
line, = ax.plot([],[],lw=2)

def init_anim():
    line.set_data([],[])
    return line,

# draw a circle with constant velocity
def animate(nitr):
    alpha = 2.*np.pi*nitr/100.
    angles = np.linspace(0,alpha,nitr+1)
    x = np.cos(angles)
    y = np.sin(angles)
    line.set_data(x,y)
    return line,

# draw a circle with constant velocity
def animate_2(nitr,x,y):
    line.set_data(x[:nitr+1],y[:nitr+1])
    return line,

# draw a circle with non-constant velocity
def animate_3(nitr):
    line.set_data(x[:nitr+1],y[:nitr+1])
    return line,

#whole_circle = np.linspace(0,2*np.pi,100)
#x = np.cos(whole_circle)
#y = np.sin(whole_circle)

n = np.arange(250)
x = np.cos(1./2.*(2.*np.pi*n/250.)**2)
y = np.sin(1./2.*(2.*np.pi*n/250.)**2)

#ani = anim.FuncAnimation(fig,animate_2,frames=len(x),init_func=init_anim,fargs=(x,y),interval=13,blit=True,repeat=False)
ani = anim.FuncAnimation(fig,animate_3,frames=len(x),init_func=init_anim,interval=13,blit=True,repeat=False)
plt.show()

plt.rcParams['animation.ffmpeg_path'] = '/usr/bin/ffmpeg'
FFwriter = anim.FFMpegWriter(fps=1000./13.)
ani.save('circle.mp4',writer=FFwriter,extra_args=['-vcodec','libx264'])
