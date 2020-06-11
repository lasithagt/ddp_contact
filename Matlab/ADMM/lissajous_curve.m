function [x,y,z] = lissajous_curve(t, z0)

%     t = linspace(-pi,pi,500);
    A = 0.1; B = 0.1; a = 1; b = 1; d = pi/2;

    x = A * sin(a*t + d);
    y = B * sin(b*t);
    z = z0 * ones(1,numel(t));
    figure(2)
    plot(x,y,'LineWidth',2)

end