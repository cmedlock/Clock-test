import math
import numpy as np
import matplotlib.pyplot as plt

n = np.linspace(0,20,100)
f0 = 1/math.sqrt(2*math.pi*3.5)*np.exp(-(n-5)**2/(2*3.5))
f1 = 1/math.sqrt(2*math.pi*3.5)*np.exp(-(n-15)**2/(2*3.5))
plt.close('all')
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(f0,color='blue',linewidth=3)
ax.plot(f1,color='red',linewidth=3)
fig.savefig('/Users/cmedlock/Documents/DSP_UROP/slides/f0_and_f1.png')

x = np.linspace(0,1,100)
y = np.sqrt(x)
ax.clear()
ax.plot(x,y,color='black',linewidth=3)
ax.fill_between(x,y,color='black',alpha=0.5)
fig.savefig('/Users/cmedlock/Documents/DSP_UROP/slides/sample_ROC.png')
plt.show()