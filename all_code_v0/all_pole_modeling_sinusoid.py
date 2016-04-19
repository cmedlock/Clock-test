import math
import numpy as np
import matplotlib.pyplot as plt

pi = math.pi
n = np.arange(250)
omega = 2.*pi/250.
x = 63.5*np.cos(omega*n)+63.5
y = 63.5*np.cos(omega*n)
shift = 15
noise = np.random.normal(loc=0,scale=2,size=len(n))
x = np.concatenate((x[shift:],x[:shift]))+noise
y = np.concatenate((y[shift:],y[:shift]))+noise

# form model using Y-W eqns

p = 3 # model order
# 'circular' autocorrelation
rxx,ryy = [],[]
x_periodic = np.concatenate((x,x))
y_periodic = np.concatenate((y,y))
for w in range(p+1):
    rxx.append(np.dot(x,x_periodic[w:w+len(x)]))
    ryy.append(np.dot(y,y_periodic[w:w+len(y)]))
    
# calculate linear prediction coefficients
D_x = np.array(rxx[1:p+1])
D_y = np.array(ryy[1:p+1])
W_x = np.empty((p,p))
W_y = np.empty((p,p))
ak_x = np.empty((p))
ak_y = np.empty((p))
#ck = []
for row in range(p):
    for column in range(row,p):
        W_x[row][column] = rxx[column-row]
        W_x[column][row] = rxx[column-row]
        W_y[row][column] = ryy[column-row]
        W_y[column][row] = ryy[column-row]
# LPC spectrum
W_x_inv = np.linalg.inv(W_x)
W_y_inv = np.linalg.inv(W_y)
#print np.dot(W_inv,W)
ak_x = np.dot(W_x_inv,D_x)
ak_y = np.dot(W_y_inv,D_y)
print sum(ak_x),sum(ak_y)
## LPC cepstrum
#ck = [ak[0]]
#for k in range(2,p+1):
#    x1 = ak[k-1]
#    for m in range(1,k):
#        x1 += float(m)/float(k)*ak[m-1]*ck[k-m-1]
#    ck.append(x1)

# impulse response of linear prediction filter
a_x,a_y = [1],[1]
for w in range(p):
    a_x.append(-float(ak_x[w]))
    a_y.append(-float(ak_y[w]))

# linear prediction error (LPE):
g_x_x = np.empty((len(x)+p)) # predict x[n] with x[n] model
g_y_y = np.empty((len(y)+p)) # predict y[n] with y[n] model
for w in range(len(x)+p):
    lpe_x_x,lpe_y_y = 0,0
    #if w<3:
        #print 'w = ',w
    for d in range(len(a_x)):
        if w-d<len(x)+p:
            #if w<3:
                #print '***   ',w-d,d,' ---> ',x_periodic[w-d],a_x[d]
            lpe_x_x += x_periodic[w-d]*a_x[d]
            lpe_y_y += y_periodic[w-d]*a_y[d]
    g_x_x[w] = lpe_x_x
    g_y_y[w] = lpe_y_y
g_x_x,g_y_y = np.array(g_x_x),np.array(g_y_y)
# linear prediction of x[n] or y[n]
xhat_x_x = -g_x_x[:-p]+x
yhat_y_y = -g_y_y[:-p]+y

plt.close('all')
fig = plt.figure()
ax1,ax2 = fig.add_subplot(211),fig.add_subplot(212)
ax1.plot(x,label='x[n]')
ax1.plot(y,label='y[n]')
ax1.legend(loc='best')
ax2.plot(g_x_x[:len(x)],label='g_x[n]')
ax2.plot(g_y_y[:len(y)],label='g_y[n]')
ax2.legend(loc='best')
#ax2.text(5,0.003,'p = '+str(p),fontsize=15)
#ax2.text(60,6*10**-10,'p = '+str(p),fontsize=15)

plt.show()
