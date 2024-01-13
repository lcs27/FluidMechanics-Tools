'''
Author: lcs27
Date: Jan. 13, 2024

This file serves to resolve the following equation with FEM
- K \frac{\partial^2 u}{\partial x^2} + r(x) = 0, x \in ]0,L[
- u(0) = u_0
- K \frac{\partial u}{\partial x}(L) = q_N

The whole system can be transformed to a linear system and we will solve this linear system

The weak formulation is 
w = u+u0
- int_0^L k w'(x) v'(x) dx + [k w'(x)v(x)]^L_0 + int_0^L r(x) v(x) dx = 0

The matrix and vector in this equations are:
- A_ij=  int_0^L k phi_i'(x) * phi_j'(x) dx 
- b_j =  int_0^L r(x) phi_j(x) dx + k*q_N phi_j(L)
'''
import numpy as np
import matplotlib.pyplot as plt

## Given of the problem
k=1
L=1
u0=1
qN=1
n=50
def r(x):
    #return 1
    #return np.sin(np.pi*x/L)
    return np.sin(np.pi*x/(2*L))

# Solution of the problem
## We use Linear basis
delta_x=L/n
K = np.zeros((n,n))
F = np.zeros((n))

## Elemental integration over each interval
def fe(x0,x1):
    # Numerical integration using trapezoidal rule
    Fe = np.array([r(x0)/2,r(x1)/2])*(x1-x0)
    return Fe
ke = np.matrix([[n/L**2,-n/L**2],[-n/L**2,n/L**2]])*k

## Assembly!
K[0,0] += ke[1,1] #Dirichlet boundary condition
for i in range(1,n):
    K[(i-1):(i+1),(i-1):(i+1)] += ke
    F[(i-1):(i+1)] += fe(i/n*L,(i+1)/n*L)

## Neumann boundary condition
F[n-1] += qN*k

## Solve linear equation
u = np.linalg.solve(K,F)+u0

# Theoretical solutions
x=np.linspace(0,L,n+1,endpoint=True)
#u_analytical = -0.5*(x**2)+2*x+u0
#u_analytical = np.sin(np.pi*x)/(np.pi**2)+(1+1/np.pi)*x+u0
u_analytical = 4*np.sin(np.pi*x/2)/(np.pi**2)+x+u0

# Plot the result
fig,ax = plt.subplots()
ax.plot(x[1:(n+1)],u,label='FEM',marker='+')
ax.plot(x,u_analytical,label='theoretical',marker='.')
ax.legend(loc='best')
ax.grid(zorder=0)
fig.savefig('FEM_Diffusion_result.png')
plt.close()