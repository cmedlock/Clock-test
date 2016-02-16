import math
import numpy as np
import matplotlib.pyplot as plt

pi = math.pi
n = np.arange(10)
omega = 2.*pi/5.
h = np.cos(omega*n)
shift = 3
h = np.concatenate((h[shift:],h[:shift]))
g = np.empty((len(h)))
a1 = 2*math.cos(omega)
a2 = -1

for w in range(len(h)):
    lin_model_error = h[w]-a1*h[w-1]-a2*h[w-2]
    g[w] = lin_model_error if lin_model_error>10**-15 else 0
    print 'lin_model_error is now ',lin_model_error

plt.close('all')
fig = plt.figure()
ax1,ax2 = fig.add_subplot(211),fig.add_subplot(212)
ax1.plot(h)
ax2.plot(g)
plt.show()
