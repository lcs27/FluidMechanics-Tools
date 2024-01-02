import numpy as np
import matplotlib.pyplot as plt
# ***********************************************************************#
# Lattice Boltzmann Method (LBM)                                         #
# original code by Li Zhe as a tutorial template and finished by LCS     #
# Email: zhe.li@ec-nantes.fr                                             #
# Version: 0.0 Dec.26 2023                                               #
# Lattice: D2Q9                                                          #
#                         6---2---5                                      #
#                         | \ | / |                                      #
#                         3---0---1                                      #
#                         | / | \ |                                      #
#                         7---4---8                                      #
#                                                                        #
# ***********************************************************************#
#
# =======================================================================#
#                               Constants                                #
# =======================================================================#
deltax = 1
deltat = 1
cs_lat = 1.0/(3**0.5)*deltax/deltat
cs2 = cs_lat**2
cs4 = cs2**2
inv1cs2 = 1.0/cs2
inv2cs2 = 1.0/(2.0*cs2)
inv1cs4 = 1.0/cs4
inv2cs4 = 1.0/(2.0*cs4)
ex = [0.0, 1.0, 0.0, -1.0, 0.0, 1.0, -1.0, -1.0, 1.0]
ey = [0.0, 0.0, 1.0, 0.0, -1.0, 1.0, 1.0, -1.0, -1.0]
w = [4.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0]
nbDist = 9
g = 5

# =======================================================================#
#                       Computational parameters                         #
# =======================================================================#
# to be modified ...
Nx = 3
Ny = 25


# =======================================================================#
#                           Initialization                               #
# =======================================================================#
fdist = np.zeros((nbDist, Nx, Ny))
fequi = np.zeros((nbDist, Nx, Ny))
Force = np.zeros((nbDist, Nx, Ny))

rho   = np.zeros((Nx, Ny))
vx    = np.zeros((Nx, Ny))
vy    = np.zeros((Nx, Ny))

rho = np.ones((Nx,Ny))
for kk in range(nbDist):
    fequi[kk,:,:] = rho*w[kk]*(1+(ex[kk]*vx+ey[kk]*vy)/cs2+(ex[kk]*vx+ey[kk]*vy)**2/2/cs4-(vx**2+vy**2)/2/cs2)

    Force[kk,:,:] = rho*w[kk]*((ex[kk]-vx)*g/cs2 + (ex[kk]*vx+ey[kk]*vy)*ex[kk]*g/cs4)

fdist=fequi
# =======================================================================#
#                             Iterations                                 #
# =======================================================================#
nbSteps = 2000
tau = 1
nu = cs2*(tau-1/2)
print('nu=',nu)
for iStep in range(0,nbSteps):
    # ---------------------------------
    # (1) Collision
    fcoll = fdist - deltat/tau*(fdist-fequi) + deltat*(1-deltat/2/tau)*Force
    
   # ---------------------------------
   # (2) Streaming
    for j in range(0,Ny):
        for i in range(0,Nx):
           for kk in range(nbDist):
                if(((j+ey[kk]) < Ny) & ((j+ey[kk]) >= 0)):
                    fdist[kk,int(np.mod(i+ex[kk],Nx)),int(j+ey[kk])] = fcoll[kk,i,j]

    
   # ---------------------------------
   # (3) Boundary condition with bounce back
   # (3.3) Bottom
    j = 0
    fdist[8,:,j] = fcoll[6,:,j]
    fdist[4,:,j] = fcoll[2,:,j]
    fdist[7,:,j] = fcoll[5,:,j]

   # (3.4) Top
    j = Ny-1
    fdist[6,:,j] = fcoll[8,:,j]
    fdist[2,:,j] = fcoll[4,:,j]
    fdist[5,:,j] = fcoll[7,:,j]
    
   # ---------------------------------
   # (4) Macroscopic variables
    rho = np.zeros((Nx, Ny))
    vx = np.zeros((Nx, Ny))
    vy = np.zeros((Nx, Ny))
    for kk in range(nbDist):
        rho = rho + fdist[kk,:,:]
        vx = vx + fdist[kk,:,:]*ex[kk] 
        vy = vy + fdist[kk,:,:]*ey[kk]
    vx = vx/rho
    vy = vy/rho
    vx = vx + deltat/2*g


    # (5) Prepare to next 
    for kk in range(nbDist):
        fequi[kk,:,:] = rho*w[kk]*(1+(ex[kk]*vx+ey[kk]*vy)/cs2+(ex[kk]*vx+ey[kk]*vy)**2/2/cs4-(vx**2+vy**2)/2/cs2)
        Force[kk,:,:] = rho*w[kk]*((ex[kk]-vx)*g/cs2 + (ex[kk]*vx+ey[kk]*vy)*ex[kk]*g/cs4)


#xx,yy = np.meshgrid(np.linspace(0,Nx,Nx), np.linspace(0,Ny,Ny))
        
# Visualisation
fig, ax = plt.subplots()
cs = ax.contourf(rho)
fig.colorbar(cs)
fig.savefig('./NumericalMethod/LBM_Poiseuille_rho.jpg',dpi=300)
plt.close()

fig, ax = plt.subplots()
cs =ax.contourf(vx)
fig.colorbar(cs)
fig.savefig('./NumericalMethod/LBM_Poiseuille_vx.jpg',dpi=300)
plt.close()



fig, ax = plt.subplots()
cs =ax.contourf(vy)
fig.colorbar(cs)
fig.savefig('./NumericalMethod/LBM_Poiseuille_vy.jpg',dpi=300)


H=(Ny-1)*deltax
y = np.linspace(0,H,Ny)
fig, ax = plt.subplots()
ax.plot(rho[0,:],label='rho')
ax.plot(vx[0,:],label='ux')
ax.plot(g*y*(H-y)/(2*nu),label='ux_ref')
ax.plot(vy[0,:],label='vy')
ax.legend()
fig.savefig('./NumericalMethod/LBM_Poiseuille_vx_profile.jpg',dpi=300)
plt.close()

