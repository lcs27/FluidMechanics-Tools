parameter_vague;
x0=0:0.25:L;
[x,y]=meshgrid(x0,-0.5:0.05:-0.05);
x1=2;
x2=2;
x3=2;
y1=-0.1;
y2=-0.2;
y3=-0.3;
for t=0:deltat:totalt
    eta=A*cos(k.*x0-omega.*t);
    area(x0,eta+y_shift,'FaceColor','c');
    hold on;
    
   [u,v]=speed_vagues(x,y,t);
    quiver(x,y+y_shift,u,v,'AutoScaleFactor',0.4,'Color','b')
    
    
    
    [u1,v1]=speed_vagues(x1,y1,t);
    x1=x1+u1*deltat;
    y1=y1+v1*deltat;
    [u2,v2]=speed_vagues(x2,y2,t);
    x2=x2+u2*deltat;
    y2=y2+v2*deltat;
    [u3,v3]=speed_vagues(x3,y3,t);
    x3=x3+u3*deltat;
    y3=y3+v3*deltat;
    scatter(x1,y1+y_shift,'r','filled');
    scatter(x2,y2+y_shift,'b','filled');
    scatter(x3,y3+y_shift,'k','filled');
    hold off;
    xlim([0,L])
    ylim([y_shift-0.5,y_shift+0.3])
    pause(0.01);
end
