function [x,y] = draw_trajectory(speedfunction,x0,y0,totalt,deltat)
x=[x0];
y=[y0];
currentx=x0;
currenty=y0;
for i=0:deltat:totalt
    [u,v] = speedfunction(currentx,currenty,i);
    currentx=currentx+deltat*u;
    currenty=currenty+deltat*v;
    x=[x,currentx];
    y=[y,currenty];
end
end

