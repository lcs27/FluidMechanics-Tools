function [x,y] = draw_trajectory(speedfunction,particle,totalt,deltat)
x=[particle.x];
y=[particle.y];
for i=0:deltat:totalt
    [u,v] = speedfunction(t,particle);
    particle.x=particle.x+deltat*u;
    particle.y=particle.y+deltat*v;
    x=[x,particle.x];
    y=[y,particle.y];
end
end

