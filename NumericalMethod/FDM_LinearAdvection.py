'''
TODO!
Author: lcs27
Date: Nov. 18, 2023

This file shows how we solve 1D linear advection equation with finite difference method.

The linearised equation is 
u_j^{n+1} = (1-a*dt/dx)*u_j^n+a*dt/dx*u_{j-1}^n
'''
import numpy as np
import matplotlib.pyplot as plt

## Configuration of the problem
a = 2
N = 256
dx = 1/N
dt = 0.001
T = 0.12
def u0(x):
    return np.exp(np.sin(2*np.pi*x)+1)
    #return np.heaviside(x-0.5,1)
    #return np.sinc(x-0.5)
x = np.linspace(0,1,N,endpoint=False)
u0s = u0(x)
u=u0s
cfl = a*dt/dx
print('CFL=',cfl)

## Solution of the problem
def update(x):
    sol = (1-a*dt/dx)*x + a*dt/dx*np.roll(x,1)
    return sol

# Go iteration!
for i in range(int(T/dt)):
    u = update(u)

## Plot the result
fig,ax = plt.subplots()
ax.plot(x,u0s,label=r'$u_0$',linestyle='--')
ax.plot(x,u,label=r'$u$')
ax.legend(loc='best')
ax.grid(zorder=0)
fig.savefig('./FDM_LinearAdvection.png')
plt.close()
