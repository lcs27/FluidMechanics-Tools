import numpy as np

# Velocity part
def advect(u,dx,dt):
    w1 = np.zeros_like(u)
    _,m,n = u.shape
    for i in range(m):
        for j in range(n):
            ii = (i + u[0,i,j]*dt/dx)%m
            jj = (j + u[1,i,j]*dt/dx)%n
            s,ii = np.modf(ii)
            t,jj = np.modf(jj)
            ii = int(ii)
            jj = int(jj)
            w1[0,ii,jj] += (1-s)*(1-t)*u[0,i,j]
            w1[0,(ii+1)%m,jj] += s*(1-t)*u[0,i,j] 
            w1[0,ii,(jj+1)%n] += (1-s)*t*u[0,i,j] 
            w1[0,(ii+1)%m,(jj+1)%n] += s*t*u[0,i,j]
            w1[1,ii,jj] += (1-s)*(1-t)*u[1,i,j]
            w1[1,(ii+1)%m,jj] += s*(1-t)*u[1,i,j] 
            w1[1,ii,(jj+1)%n] += (1-s)*t*u[1,i,j] 
            w1[1,(ii+1)%m,(jj+1)%n] += s*t*u[1,i,j]
    print(w1[0,0,0])
    return w1

def advectbis(u,dx,dt):
    w1 = np.copy(u)
    _,m,n = u.shape
    for _ in range(20):
        for i in range(m):
            for j in range(n):
                ii = (i - w1[0,i,j]*dt/dx)%m
                jj = (j - w1[1,i,j]*dt/dx)%n
                s,ii = np.modf(ii)
                t,jj = np.modf(jj)
                ii = int(ii)
                jj = int(jj)
                w1[0,i,j] = (1-s)*(1-t)*u[0,ii,jj] + s*(1-t)*u[0,(ii+1)%m,jj] + (1-s)*t*u[0,ii,(jj+1)%n] + s*t*u[0,(ii+1)%m,(jj+1)%n]
                w1[1,i,j] = (1-s)*(1-t)*u[1,ii,jj] + s*(1-t)*u[1,(ii+1)%m,jj] + (1-s)*t*u[1,ii,(jj+1)%n] + s*t*u[1,(ii+1)%m,(jj+1)%n]
    return w1

def diffuse(w1,u,Re,dt,dx):
    w2 = np.zeros_like(w1)
    _,m,n = u.shape
    for i in range(m):
        for j in range(n):
            w2[0,i,j] = w1[0,i,j] + dt/(Re*dx*dx) * (-4*u[0,i,j]+u[0,(i-1)%m,j]+u[0,(i+1)%m,j]+u[0,i,(j-1)%n]+u[0,i,(j+1)%n])
            w2[1,i,j] = w1[1,i,j] + dt/(Re*dx*dx) * (-4*u[1,i,j]+u[1,(i-1)%m,j]+u[1,(i+1)%m,j]+u[1,i,(j-1)%n]+u[1,i,(j+1)%n])
    print('w1',w1[0,0,0],'w2',w2[0,0,0])
    return w2

def diffusebis(w1,Re,dt,dx):
    w2 = np.copy(w1)
    _,m,n = w1.shape
    for _ in range(20):
        for i in range(m):
            for j in range(n):
                w2[0,i,j] = (w1[0,i,j] + dt/(Re*dx*dx) * (w2[0,(i-1)%m,j]+w2[0,(i+1)%m,j]+w2[0,i,(j-1)%n]+w2[0,i,(j+1)%n]))/(1+4*dt/(Re*dx*dx))
                w2[1,i,j] = (w1[1,i,j] + dt/(Re*dx*dx) * (w2[1,(i-1)%m,j]+w2[1,(i+1)%m,j]+w2[1,i,(j-1)%n]+w2[1,i,(j+1)%n]))/(1+4*dt/(Re*dx*dx))
    return w2

def addExternalForces(w2,f,dt):
    return w2+f*dt

def project(w3,dx):
    _,m,n = w3.shape
    p = np.zeros((m,n))
    theta = np.zeros((m,n))
    for i in range(m):
        for j in range(n):
            theta[i,j] = 0.5 *(w3[0,(i+1)%m,j]-w3[0,(i-1)%m,j] + w3[1,i,(j+1)%n]-w3[1,i,(j-1)%n]) / dx
    for _ in range(20):
        for i in range(m):
            for j in range(n):
                p[i,j] = 0.25*(- theta[i,j] * (dx**2) + p[(i-1)%m,j] + p[(i+1)%m,j] + p[i,(j-1)%n] + p[i,(j+1)%n]) 
    for i in range(m):
        for j in range(n):
            w3[0,i,j] -= 0.5 * (p[(i+1)%m,j]-p[(i-1)%m,j]) / dx
            w3[1,i,j] -= 0.5 * (p[i,(j+1)%n]-p[i,(j-1)%n]) / dx
    return w3,p


def compute_energy(u,dx):
    E = 0
    d,m,n = u.shape
    for k in range(d):
        for i in range(m):
            for j in range(n):
                E += 1/2*u[k,i,j]**2 * dx*dx
    return E


# Passive scalar part
def advect_scalar(a,u,dt,dx):
    a1 = np.zeros_like(a)
    m,n = a.shape
    for i in range(m):
        for j in range(n):
            ii = (i + u[0,i,j]*dt/dx)%m
            jj = (j + u[1,i,j]*dt/dx)%n
            s,ii = np.modf(ii)
            t,jj = np.modf(jj)
            ii = int(ii)
            jj = int(jj)
            a1[ii,jj] += (1-s)*(1-t)*a[i,j]
            a1[(ii+1)%m,jj] += s*(1-t)*a[i,j] 
            a1[ii,(jj+1)%n] += (1-s)*t*a[i,j] 
            a1[(ii+1)%m,(jj+1)%n] += s*t*a[i,j]
    return a1

def diffuse_scalar(a1,a,kappa,dt,dx):
    a2 = np.zeros_like(a1)
    m,n = a1.shape
    for i in range(m):
        for j in range(n):
            a2[i,j] = a1[i,j] + kappa*dt/(dx*dx) * (-4*a[i,j]+a[(i-1)%m,j]+a[(i+1)%m,j]+a[i,(j-1)%n]+a[i,(j+1)%n])
    return a2

def dissipate_scalar(a2,a,beta,dt):
    return a2-beta*a*dt

def source_scalar(a3,S,dt):
    return a3+S*dt


