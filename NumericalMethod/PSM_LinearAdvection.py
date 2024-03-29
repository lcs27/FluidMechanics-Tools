'''
TODO!
Author: lcs27
Date: Nov. 18, 2023

This file shows how we solve 1D linear advection equation with spectral method.

The linearised equation is 
uhat^{n+1}_k=(1-aik\delta t)uhat^n_k
'''
import numpy as np
import matplotlib.pyplot as plt

## Configuration of the problem
a = 2
N = 256
dx = 1/N
dt = 0.0002
T = 0.12
def u0(x):
    #return np.log(np.abs(x-0.5)+1)
    #return np.heaviside(x-0.5,1)
    return np.exp(np.sin(2*np.pi*x)+1)
x = np.linspace(0,1,N,endpoint=False)
u0s = u0(x)
u=u0s
ks = np.fft.fftfreq(N,dx)*(2*np.pi)
cfl = a*dt*np.max(ks)
print('CFL=',cfl)

## Solution of the problem
def update(x):
    xf = np.fft.fft(x)/np.sqrt(N)
    sol = (1-a*1j*ks*dt)*xf
    return np.real(np.fft.ifft(sol)*np.sqrt(N))

# Go iteration!
for i in range(int(T/dt)):
    u = update(u)

## Plot the result
fig,ax = plt.subplots()
ax.plot(x,u0s,label=r'$u_0$',linestyle='--')
ax.plot(x,u,label=r'$u$')
ax.legend(loc='best')
ax.grid(zorder=0)
fig.savefig('./LinearAdvectionPSM.png')
plt.close()




