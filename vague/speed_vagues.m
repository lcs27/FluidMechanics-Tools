function [u,v] = speed_vagues(x,y,t)
parameter_vague;
u=A*omega*cosh(k.*y+k*h)./sinh(k*h).*cos(k.*x-omega.*t);
v=A*omega*sinh(k.*y+k*h)./sinh(k*h).*sin(k.*x-omega.*t);
end

