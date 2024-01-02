import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
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

## Physical Parameters
deltax_phy = 0.1
rho_phy = 1
g_phy = 10
nu_phy = 4e-1
Ly_phy = 1
nbSteps = 1200
BC = 1 # 0 = Zou-He, 1= Bounce back

## Rescaling parameters
deltat_phy = deltax_phy**2/nu_phy/6
deltax = 1
deltat = 1
Cx = deltax_phy/deltax
Ct = deltat_phy/deltat
Crho = 1
Cu = Cx/Ct
Cp = Crho*Cu**2
Cg = Cx/(Ct**2)
Cnu = Cx**2/Ct

cs = 1.0/(3**0.5)*deltax/deltat
cs2 = cs**2
cs4 = cs2**2
g = g_phy/Cg
nu = nu_phy/Cnu
tau = (nu/cs2+1/2)

ex = [0.0, 1.0, 0.0, -1.0, 0.0, 1.0, -1.0, -1.0, 1.0]
ey = [0.0, 0.0, 1.0, 0.0, -1.0, 1.0, 1.0, -1.0, -1.0]
w = [4.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0]
nbDist = 9

# =======================================================================#
#                       Computational parameters                         #
# =======================================================================#
Nx = 3
Ny = int(Ly_phy/Cx)+1
Ly_phy = (Ny-1)*Cx*deltax
print('Ly=',Ly_phy)


# =======================================================================#
#                           Initialization                               #
# =======================================================================#
y_phy = np.linspace(0,Ly_phy,Ny)
fdist = np.zeros((nbDist, Nx, Ny))
fequi = np.zeros((nbDist, Nx, Ny))
Force = np.zeros((nbDist, Nx, Ny))

rho   = np.zeros((Nx, Ny))
vx    = np.zeros((Nx, Ny))
vy    = np.zeros((Nx, Ny))

rho = np.ones((Nx,Ny))*rho_phy/Crho
for kk in range(nbDist):
    fequi[kk,:,:] = rho*w[kk]*(1+(ex[kk]*vx+ey[kk]*vy)/cs2+(ex[kk]*vx+ey[kk]*vy)**2/2/cs4-(vx**2+vy**2)/2/cs2)

    Force[kk,:,:] = rho*w[kk]*((ex[kk]-vx)*g/cs2 + (ex[kk]*vx+ey[kk]*vy)*ex[kk]*g/cs4)

fdist=fequi
# =======================================================================#
#                             Iterations                                 #
# =======================================================================#
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
    if(BC==0):
        # (3) No slip boundary condition with bounce back
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
    
    elif(BC==1):
        # (3bis) No slip boundary condition with Zou-he
        # (3bis.3) Bottom
        j = 0
        fdist[4,:,j] = fcoll[2,:,j]
        fdist[7,:,j] = fcoll[5,:,j] + 0.5*(fcoll[1,:,j]-fcoll[3,:,j])
        fdist[8,:,j] = fcoll[6,:,j] - 0.5*(fcoll[1,:,j]-fcoll[3,:,j])

        # (3bis.4) Top
        j = Ny-1
        fdist[2,:,j] = fcoll[4,:,j]
        fdist[5,:,j] = fcoll[7,:,j] - 0.5*(fcoll[1,:,j]-fcoll[3,:,j])
        fdist[6,:,j] = fcoll[8,:,j] + 0.5*(fcoll[1,:,j]-fcoll[3,:,j])
        
    else:
        print("error")


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


# print('Cu=',Cu)
fig, ax = plt.subplots()
ax.plot(y_phy,rho[0,:]*Crho,label='rho')
ax.plot(y_phy,vx[0,:]*Cu,label='ux',marker='+')
ax.plot(y_phy,g_phy*y_phy*(Ly_phy-y_phy)/(2*nu_phy),label='ux_ref',marker='x')
ax.plot(y_phy,vy[0,:],label='vy')
ax.set_title('t='+str(deltat_phy*nbSteps))
ax.legend()
fig.savefig('./NumericalMethod/LBM_Poiseuille_vx_profile.jpg',dpi=300)
plt.close()

