clc;
a = 20;
b = 8/3;
c = 28;
x0 = 1;
y0 = 5;
z0 = 9;
time = [0, 50000];


[t, XYZ] = ode45(@(t, X) lorentz(t, X, a, b, c), time, [x0; y0; z0]);

x = XYZ(:, 1);
y = XYZ(:, 2);
z = XYZ(:, 3);
xy = x .*y ;

cum_integral = cumtrapz(t, xy);       
durations = t - t(1);               


running_avg = zeros(size(t));
running_avg(1) = xy(1);              
running_avg(2:end) = cum_integral(2:end) ./ durations(2:end);


figure;
plot(t, running_avg, 'r', 'LineWidth', 1.5);
xlabel('Time');
ylabel('Average Value of x');
title('Running Time-Average of x in Lorenz System');
grid on;

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