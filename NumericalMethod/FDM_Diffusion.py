'''
Author: lcs27
Date: Nov. 5, 2023

This file serves to resolve the following equation with FDM
- K \frac{\partial^2 u}{\partial x^2} + r(x) = 0, x \in ]0,L[
- u(0) = 0
- K \frac{\partial u}{\partial x}(L) = q_N

The whole system can be transformed to a linear system and we will solve this linear system, we will then compare the result of backward and centred method 
'''
import numpy as np
import matplotlib.pyplot as plt

## Given of the problem
K=1
L=1
N=50
delta_x=L/N
u0=0
qN=1
def r(x):
    #return 1
    #return np.sin(np.pi*x/L)
    return np.sin(np.pi*x/(2*L))

## Stage 1
A = np.zeros((N+1,N+1))
b = np.zeros((N+1,1))
for i in range(1,N):
    A[i,i-1]=1
    A[i,i]=-2
    A[i,i+1]=1
for i in range(1,N):
    b[i]=-r(i*delta_x)*(delta_x**2)/K

## Stage 2
A[0,0]=1
b[0]= u0
A[N,N-1]=-1
A[N,N]=1
b[N]=qN*delta_x/K
u = np.squeeze(np.linalg.solve(A, b))

print(A)

## Stage 2-bis, use fictitious point method
Abis = np.zeros((N+2,N+2))
bbis = np.zeros((N+2,1))

# A fictitious point of x_{N+1}
for i in range(1,N+1):
    Abis[i,i-1]=1
    Abis[i,i]=-2
    Abis[i,i+1]=1
for i in range(1,N+1):
    bbis[i]=-r(i*delta_x)*(delta_x**2)/K
Abis[0,0]=1
bbis[0]= u0
Abis[N+1,N-1]=-1
Abis[N+1,N+1]=1
bbis[N+1]=qN/K*2*delta_x

print(Abis)

ubis = np.squeeze(np.linalg.solve(Abis, bbis))

## Analytical solution
x=np.array(range(0,N+1))/N*L
# u_analytical = -0.5*(x**2)+2*x
#u_analytical = np.sin(np.pi*x)/(np.pi**2)+(1+1/np.pi)*x
u_analytical = 4*np.sin(np.pi*x/2)/(np.pi**2)+x

## Plot the result
fig, ax = plt.subplots()
ax.plot(x,u,label='backward',marker='+')
ax.plot(x,ubis[range(0,N+1)],label='centered',marker='x')
ax.plot(x,u_analytical,label='theoretical',marker='.')
ax.legend()
fig.savefig('FDM_Diffusion_result.png', dpi=300)
plt.close()

## Calculate the error
error = max(np.abs(u_analytical-u))
error_bis = max(np.abs(u_analytical-ubis[range(0,N+1)]))

print(error,error_bis)
fig, ax = plt.subplots()
ax.semilogy(x,np.abs(u_analytical-u),label='backward',marker='x')
ax.semilogy(x,np.abs(u_analytical-ubis[range(0,N+1)]),label='centered',marker='+')
ax.legend()
fig.savefig('FDM_Diffusion_error.png', dpi=300)
plt.close()