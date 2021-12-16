function [u,v] = speed_interfaces(t,particle)
param_interface;
x=particle.x;
y=particle.y;

if particle.type == 0
    u=U1+real(1i*k*alpha1*e^(k*y)*e^(1i*w*t-1i*k*x));
    v=real(k*alpha1*e^(k*y)*e^(1i*w*t-1i*k*x));
else
    u=U2+real(1i*k*beta2*e^(-k*y)*e^(1i*w*t-1i*k*x));
    v=real(-k*beta2*e^(-k*y)*e^(1i*w*t-1i*k*x));
end

