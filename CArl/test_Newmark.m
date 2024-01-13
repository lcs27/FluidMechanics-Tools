N=400;
K=zeros(N,N);
M=zeros(N,N);
E=200000;
rho=1000;
L=6;


for ii=1:N
    K(ii,ii)=2*E*N/L;
    M(ii,ii)=2*rho*L./(3*N);

    if (ii-1)>0
        K(ii-1,ii)=-E*N/L;
        M(ii-1,ii)=rho*L./(6*N);
    end

    if (ii+1)<=N
        K(ii+1,ii)=-E*N/L;
        M(ii+1,ii)=rho*L./(6*N);
    end
end

F=zeros(N,1);
F(N)=1;

beta=1/4;
gamma=1/2;
dt=0.01;
s_t=10*dt;

A=M+beta*dt^2*K;

Udd=zeros(N,1);
Ud=zeros(N,1);
U=zeros(N,1);

pU=zeros(N,1);
pUd=zeros(N,1);

times=s_t/dt;

us=zeros(N,times);

for kkk=1:times

    % This is moment kkk-1
    pU=U+dt*Ud+dt^2*(1/2-beta)*Udd;
    pUd=Ud+dt*(1-gamma)*Udd;

    % This is moment kkk
    B=F-K*pU;
    Udd=A\B;

    U=pU+beta*dt^2*Udd;
    Ud=pUd+gamma*dt*Udd;

    us(:,kkk)=U;
end
%plot((1:times)*dt,us)
plot(U)
