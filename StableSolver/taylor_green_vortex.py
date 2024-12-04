import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import PillowWriter
from solver import *

# Initialization
N = 40
L = 2*np.pi
dx = L/N
dt = 0.04
alpha = 0
beta = 0.01
Re = 100
T = 1
kappa = 0.01
x = np.meshgrid(np.linspace(dx/2,L-dx/2,N,True),np.linspace(dx/2,L-dx/2,N,True))
u = np.array([np.sin(x[0]) * np.cos(x[1]),- np.cos(x[0]) * np.sin(x[1])])
f = np.array([alpha*np.sin(x[0]) * np.cos(x[1]),-alpha*np.cos(x[0]) * np.sin(x[1])])
S = np.zeros((N,N))
a = np.random.rand(N,N)
E = []
t = []

# Animation preparation
metadata = dict(title='TGV', artist='Matplotlib',comment='TGV')
writer = PillowWriter(fps=15, metadata=metadata)
metadata1 = dict(title='TGV scalar', artist='Matplotlib',comment='TGV scalar')
writer1 = PillowWriter(fps=15, metadata=metadata1)
fig, ax = plt.subplots(1,2,figsize=(6,3))
levels = np.linspace(np.min(u),np.max(u),100,endpoint=True)
ax[0].contourf(x[0],x[1],u[0],levels)
ax[1].contourf(x[0],x[1],u[1],levels)
fig.suptitle(r'$t='+str(0)+r'$')
fig1, ax1 = plt.subplots()
levels1 = np.linspace(0,1,100,endpoint=True)
ax1.contourf(x[0],x[1],a,levels1)
fig1.suptitle(r'$t='+str(0)+r'$')
fig.savefig("./StableSolver/u0.jpg")
E.append(compute_energy(u,dx))
t.append(0)

# Solving loop
with writer.saving(fig, "./StableSolver/u.gif", dpi=100):
    writer.grab_frame()
    with writer1.saving(fig1, "./StableSolver/a.gif", dpi=100):
        writer1.grab_frame()
        for i in range(int(T/dt)):
            a1 = advect_scalar(a,u,dt,dx)
            a2 = diffuse_scalar(a1,a,kappa,dt,dx)
            a3 = dissipate_scalar(a2,a,beta,dt)
            a = source_scalar(a3,S,dt)
            ax1.contourf(x[0],x[1],a,levels1)
            fig1.suptitle(r'$t='+'%.2f'%((i+1)*dt)+r'$')
            writer1.grab_frame()

            w1 = advectbis(u,dx,dt)
            w2 = diffusebis(w1,Re,dt,dx)
            w3 = addExternalForces(w2,f,dt)
            u,p = project(w3,dx)
            ax[0].contourf(x[0],x[1],u[0],levels)
            ax[1].contourf(x[0],x[1],u[1],levels)
            fig.suptitle(r'$t='+'%.2f'%((i+1)*dt)+r'$')
            writer.grab_frame()
            E.append(compute_energy(u,dx))
            t.append(dt*(i+1))
fig.savefig("./StableSolver/uT.jpg")

# Plot E-t figure
fig2, ax2= plt.subplots()
ax2.semilogy(t,E)
fig2.savefig("./StableSolver/E.jpg")


