import math
import numpy as np
import matplotlib.pyplot as plt

pi = math.pi
n = np.arange(10)
omega = 2.*pi/5.
h = np.cos(omega*n)
shift = 3
h = np.concatenate((h[shift:],h[:shift]))

# form model using Y-W eqns

p = 2 # model order
# 'circular' autocorrelation
rhh = []
h_periodic = np.concatenate((h,h))
for w in range(p+1):
    rhh.append(np.dot(h,h_periodic[w:w+len(h)]))
    
# calculate linear prediction coefficients
D = np.array(rhh[1:p+1])
W = np.empty((p,p))
ak = np.empty((p))
ck = []
for row in range(p):
    for column in range(row,p):
        W[row][column] = rhh[column-row]
        W[column][row] = rhh[column-row]
        # LPC spectrum
        W_inv = np.linalg.inv(W)
        ak = np.dot(W_inv,D)
        # LPC cepstrum
        ck = [ak[0]]
        for k in range(2,p+1):
            x1 = ak[k-1]
            for m in range(1,k):
                x1 += float(m)/float(k)*ak[m-1]*ck[k-m-1]
            ck.append(x1)

# check
g = np.empty((len(h)))
for w in range(len(h)):
    lin_model_error = h[w]
    for d in range(len(ak)):
        lpc = ak[d]
        lin_model_error -= lpc*h[w-d-1]
    g[w] = lin_model_error if lin_model_error>10**-15 else 0
    print 'lin_model_error is now ',lin_model_error


plt.close('all')
fig = plt.figure()
ax1,ax2 = fig.add_subplot(211),fig.add_subplot(212)
ax1.plot(h)
ax2.plot(g)
plt.show()
