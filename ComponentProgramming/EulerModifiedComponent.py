'''
This file has not been finished
Author: lcs27
Date: Dec. 2, 2023

This file shows a component solution and a non-component solution of ODE

x'(t) = y(t) + 2exp(t) = f(x,y,t)
y'(t) = x(t) + 1 = g(x,y,t)
x(0) = 2, y(0)= 0 
The analytical solution is 
x(t) = exp(t) + exp(-t) + (t+1)exp(t) - 1
y(t) = exp(t) - exp(-t) + t*exp(t)
'''
import numpy as np
import matplotlib.pyplot as plt

def f(x,y,t):
    return y + 2*np.exp(t)

def g(x,y,t):
    return x + 1

def A(X,t):
    return np.array([X[1]+2*np.exp(t),X[0]+1])

def analytical(t):
    return np.array([1*np.exp(t) + 2*np.exp(-t) + (t+1)*np.exp(t), 1*np.exp(t) - 2*np.exp(-t) + t*np.exp(t) ])

dts = [1e-2,1e-3,5e-4,1e-4]
t = 3


## Non-component solution
# X=[x,y]
# Advance one step
def ModifiedEulerNonComponent(X,A,t,dt):
    step1= X+dt*A(X,t)
    return X + dt/2*(A(X,t)+A(step1,t+dt))

errorx = np.zeros_like(dts)
errory = np.zeros_like(dts)
for i in range(np.size(dts)):
    X = analytical(t)
    X = ModifiedEulerNonComponent(X,A,t,dts[i])
    errorx[i] = np.abs(X[0]-analytical(t+dts[i])[0])
    errory[i] = np.abs(X[1]-analytical(t+dts[i])[1])

# Plot the result
fig,ax = plt.subplots()
ax.loglog(dts,errorx,marker='+',label=r'$e_x$, Non-component')
ax.loglog(dts,errory,marker='x',label=r'$e_y$, Non-component')

## Component solution
def ModifiedEulerComponentWithCommunicate(X,f,g,t,dt):
    x0 = X[0]
    y0 = X[1]

    # Step 1
    fx1 = f(x0,y0,t)
    fy1 = g(x0,y0,t)
    x1 = x0+t*fx1
    y1 = y0+t*fy1
    
    # Step 2
    fx2 = f(x1,y1,t+dt)
    fy2 = g(x1,y1,t+dt)
    return np.array([x0+dt/2*(fx1+fx2),y0+dt/2*(fy1+fy2)])

errorx = np.zeros_like(dts)
errory = np.zeros_like(dts)
for i in range(np.size(dts)):
    X = analytical(t)
    X = ModifiedEulerComponentWithCommunicate(X,f,g,t,dts[i])
    errorx[i] = np.abs(X[0]-analytical(t+dts[i])[0])
    errory[i] = np.abs(X[1]-analytical(t+dts[i])[1])

ax.loglog(dts,errorx,label=r'$e_x$, Component with communicate')
ax.loglog(dts,errory,label=r'$e_y$, Component with communicate')


def ModifiedEulerComponentWithoutCommunicate(X,f,g,t,dt):
    x0 = X[0]
    y0 = X[1]

    # Step 1
    fx1 = f(x0,y0,t)
    fy1 = g(x0,y0,t)
    x1 = x0+t*fx1
    y1 = y0+t*fy1
    
    # Step 2
    fx2 = f(x1,y0,t+dt)
    fy2 = g(x0,y1,t+dt)
    return np.array([x0+dt/2*(fx1+fx2),y0+dt/2*(fy1+fy2)])

errorx = np.zeros_like(dts)
errory = np.zeros_like(dts)
for i in range(np.size(dts)):
    X = analytical(t)
    X = ModifiedEulerComponentWithoutCommunicate(X,f,g,t,dts[i])
    errorx[i] = np.abs(X[0]-analytical(t+dts[i])[0])
    errory[i] = np.abs(X[1]-analytical(t+dts[i])[1])

ax.loglog(dts,errorx,label=r'$e_x$, Component without communicate')
ax.loglog(dts,errory,label=r'$e_y$, Component without communicate')
ax.legend(loc='best')
ax.grid(zorder=0)
fig.savefig('./Error.png')
plt.close()