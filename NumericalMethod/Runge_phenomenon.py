'''
Author: lcs27
Date: Nov. 5, 2023

This file shows the runge's phenomenon and its remedies
'''
import numpy as np
import matplotlib.pyplot as plt

def y(x):
    return 1/(1+(5*x)**2)


x = np.arange(-1,1,0.01)

# interpolater with n order polynomial
# a0 + a1 x + a2 x^2 + ... + an x^n = y(x), x taken as -1+2*i/n, i in [[0,n]]

def polynomial_interpolate(n):
    A=np.zeros((n+1,n+1))
    b=np.zeros((n+1,))
    for i in range(n+1):
        for j in range(n+1):
            A[i,j]=(-1+2*i/n)**j

        b[i]=y(-1+2*i/n)
    return np.flipud(np.linalg.solve(A, b))

fig, ax = plt.subplots()
ax.plot(x,y(x), label='original')
ax.plot(x, np.polyval(polynomial_interpolate(5),x), label=r'$P_5(x)$')
ax.plot(x, np.polyval(polynomial_interpolate(7),x), label=r'$P_7(x)$')
ax.plot(x, np.polyval(polynomial_interpolate(9),x), label=r'$P_9(x)$')
ax.legend()
fig.savefig('runge_phenomenon.jpg', dpi=300)


## How to remedy?
# cf. https://en.wikipedia.org/wiki/Runge%27s_phenomenon#Mitigations

# Method 1: Change of interpolation points to chebyshev nodes
# Detail source: https://en.wikipedia.org/wiki/Chebyshev_nodes
def polynomial_interpolate_remedy1(n):
    A=np.zeros((n+1,n+1))
    b=np.zeros((n+1,))
    for i in range(n+1):
        if i != 0 and i != n:
            node = np.cos((2*i-1)/2/(n-1)*np.pi)
        elif i==0:
            node = -1
        else:
            node = 1
        for j in range(n+1):
            A[i,j]=(node)**j

        b[i]=y(node)
    return np.flipud(np.linalg.solve(A, b))

fig, ax = plt.subplots()
ax.plot(x,y(x), label='original')
ax.plot(x, np.polyval(polynomial_interpolate_remedy1(5),x), label=r'$P_5(x)$',linestyle = '--')
ax.plot(x, np.polyval(polynomial_interpolate_remedy1(6),x), label=r'$P_6(x)$',linestyle = '--')
ax.plot(x, np.polyval(polynomial_interpolate_remedy1(7),x), label=r'$P_7(x)$',linestyle = '--')
ax.plot(x, np.polyval(polynomial_interpolate_remedy1(8),x), label=r'$P_8(x)$',linestyle = '--')
ax.plot(x, np.polyval(polynomial_interpolate_remedy1(9),x), label=r'$P_9(x)$',linestyle = '--')
ax.legend()
fig.savefig('runge_phenomenon_remedy1.jpg', dpi=300)

# Method 2: Least squares fitting
# Detail source: https://mathworld.wolfram.com/LeastSquaresFittingPolynomial.html
def polynomial_interpolate_remedy2(n):
    # TODO
    return np.flipud(np.linalg.solve(A, b))