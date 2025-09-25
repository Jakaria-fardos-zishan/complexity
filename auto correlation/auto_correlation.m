clc;
clear;

a = 10;
b = 8/3;
c = 28;
x0 = 1;
y0 = 5;
z0 = 9;

time = [0, 10000];
[t, XYZ] = ode45(@(t, X) lorentz(t, X, a, b, c), time, [x0; y0; z0]);

x = XYZ(:, 1);
y = XYZ(:, 2);
z = XYZ(:, 3);

maxLag = 200;

r_x = autocorr(x, maxLag);
r_y = autocorr(y, maxLag);
r_z = autocorr(z, maxLag);

lags = 0:maxLag;
figure;
subplot(3,1,1);
plot(lags, r_x, 'b-', 'LineWidth', 1.5);
xlabel('Lag'); ylabel('r_x(k)');
title('Autocorrelation of x(t)');
grid on;

subplot(3,1,2);
plot(lags, r_y, 'r-', 'LineWidth', 1.5);
xlabel('Lag'); ylabel('r_y(k)');
title('Autocorrelation of y(t)');
grid on;

subplot(3,1,3);
plot(lags, r_z, 'g-', 'LineWidth', 1.5);
xlabel('Lag'); ylabel('r_z(k)');
title('Autocorrelation of z(t)');
grid on;

function dXdt = lorentz(~, X, a, b, c)
    x = X(1); y = X(2); z = X(3);
    dx = a * (y - x);
    dy = c * x - x * z - y;
    dz = x * y - b * z;
    dXdt = [dx; dy; dz];
end

function r = autocorr(series, maxLag)
    N = length(series);
    mu = mean(series);
    denom = sum((series - mu).^2);
    r = zeros(maxLag+1, 1);
    for k = 0:maxLag
        num = sum((series(1:N-k) - mu) .* (series(1+k:N) - mu));
        r(k+1) = num / denom;
    end
end
