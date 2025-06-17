clc;
a = 20;
b = 8/3;
c = 28;
x0 = 1;
y0 = 5;
z0 = 9;
time = [0, 5000];

[t, XYZ] = ode45(@(t, X) lorentz(t, X, a, b, c), time, [x0; y0; z0]);

x = XYZ(:, 1);
y = XYZ(:, 2);
z = XYZ(:, 3);
xy = x .* y;

integral_x = trapz(t, xy);       
total_time = t(end) - t(1);    
average_x = integral_x / total_time;

fprintf('Average value of x: %.6f\n', average_x);

% Lorenz system equations
function dXdt = lorentz(~, X, a, b, c)
    x = X(1);
    y = X(2);
    z = X(3);
    
    dx = a * (y - x);
    dy = c * x - x * z - y;
    dz = x * y - b * z;
    
    dXdt = [dx; dy; dz];
end