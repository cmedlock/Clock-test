import math
import matplotlib
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as anim
import numpy as np
import os
from pylab import *

def update_lines(num, dataLines, lines):
    for line, data in zip(lines, dataLines):
        # NOTE: there is no .set_data() for 3 dim data...
        line.set_data(data[0:2, :num])
        line.set_3d_properties(data[2, :num])
    return lines

plt.close('all')

# Attaching 3D axis to the figure
fig = plt.figure()
ax = p3.Axes3D(fig)

# 3-D line of a circle drawn at non-constant velocity
n = np.arange(250)
x = np.cos(1./2.*(2.*np.pi*n/250.)**2)
y = np.sin(1./2.*(2.*np.pi*n/250.)**2)
data = [np.array([x,y,n/100.])]

# Creating line object
# NOTE: Can't pass empty arrays into 3d version of plot()
lines = [ax.plot(dat[0, 0:1], dat[1, 0:1], dat[2, 0:1])[0] for dat in data]

# Setting the axes properties
ax.set_xlim3d([-1.5, 1.5])
ax.set_xlabel('x')

ax.set_ylim3d([-1.5, 1.5])
ax.set_ylabel('y')

ax.set_zlim3d([-1.0, 3.0])
ax.set_zlabel('n')

ax.set_title('3D Test')

# Creating the Animation object
ani = anim.FuncAnimation(fig, update_lines, 250, fargs=(data, lines), interval=50, blit=False)

plt.show()

plt.rcParams['animation.ffmpeg_path'] = '/usr/bin/ffmpeg'
FFwriter = anim.FFMpegWriter(fps=1000./13.)
ani.save('circle_3d.mp4',writer=FFwriter,extra_args=['-vcodec','libx264'])
