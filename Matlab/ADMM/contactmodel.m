function []= contactmodel(R,F,E,n)
%This function is used to plot the indentation plot of contact between a
%rigid sphere and an elastic half space
%R: The Radius of ball
%F: The applied load force in ball
%E: The elastic moduli of surface
%n: The Poisson's ratio of surface
close all
a=(0.75*F*R*(1-n^2)./E)^(1/3);l
pm=F./(pi*a^2);
k=1.5*(1-n^2)./E*pm;%coefficient of curvance function, no pyhsical meaning.
%f1=@(x,y) -k*pi*0.25*a*(2-x^2-y^2);
%fsurf(f1,[-1 1 -1 1]);
funx1 = @(u,v) a*u.*sin(v);%u=r/a v=theta
funy1 = @(u,v) a*u.*cos(v);
funz1 = @(u,v) -k*pi*0.25*a*(2-u^2);
fsurf(funx1,funy1,funz1,[0 1 0 2*pi]);
%f2=@(x,y) -k*0.5*a*[(2-x^2-y^2)*asin(1./(x^2+y^2)^0.5)+(x^2+y^2-1)^0.5];
%fsurf(f2,[-3 3 -3 3]);
hold on;
funx2 = @(u,v) a*u.*sin(v);
funy2 = @(u,v) a*u.*cos(v);
funz2 = @(u,v) -k*0.5*a*[(2-u^2)*asin(1./(u^2)^0.5)+(u^2-1)^0.5];
fsurf(funx2,funy2,funz2,[1 5 0 2*pi]);
end

