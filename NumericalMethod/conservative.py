'''
Author: lcs27
Date: Dec. 9, 2024

This file shows the difference between the conservative form and non-conservative form in numerical simulation.
'''
import numpy as np
import matplotlib.pyplot as plt
N = 40
dx = 1/N
dt = 0.04
u = np.random.rand(N)

# Non-conservative
u1 = np.zeros_like(u)
for i in range(N):
    u1[i] = u[i] - dt * u[i] * (u[(i+1)%N]-u[i])/dx

# Conservative
u2 = np.zeros_like(u)
for i in range(N):
    u2[i] = u[i] - dt * (u[(i+1)%N]**2 - u[i]**2)/2/dx

print(np.sum(u),np.sum(u1),np.sum(u2))