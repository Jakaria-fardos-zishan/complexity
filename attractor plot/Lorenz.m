clc;
a = 10;
b = 5;
c = 28;
x0 = 1;
y0 = 5;
z0 = 9;

time = [0, 10000];

[t, XYZ] = ode45(@(t, X) lorentz(t, X, a, b, c), time, [x0; y0; z0]);

x = XYZ(:, 1);
y = XYZ(:, 2);
z = XYZ(:, 3);

figure
plot3(x,y,z, 'k-');

xlabel('X');
ylabel('Y');
zlabel('Z');
title('Lorenz Attractor');

function dXdt = lorentz(~, X, a, b, c)
    x = X(1);
    y = X(2);
    z = X(3);
    
    dx = a*(y - x);
    dy = c*x - x*z - y;
    dz = x*y - b*z;
    
    dXdt = [dx; dy; dz];
end