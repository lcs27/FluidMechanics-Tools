rho1=10^3;
rho2=1.1*10^3;
U1=2;
U2=3;

p0=10^5;
g=9.8;

c=1;
k=1;
w=k*(rho1*U1+rho2*U2)/(rho1+rho2)+sqrt(-k^2*rho1*rho2*(U1-U2)^2/(rho1+rho2)^2+k*g*(rho1-rho2)/(rho1+rho2));

etaA=0.1;
alpha1=1i*etaA*(U1-c);
beta2=-1i*etaA*(U2-c);

