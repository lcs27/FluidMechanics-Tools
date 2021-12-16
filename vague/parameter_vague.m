h=1;
g=9.8;
c=0.98*sqrt(g*h);
syms f;
k=0.3532
%k=double(vpasolve(f*c^2/g-tanh(f*h),f,20));
omega=k*c;
A=0.05;
L=10;
deltat=0.1;
totalt=10;
y_shift=5;