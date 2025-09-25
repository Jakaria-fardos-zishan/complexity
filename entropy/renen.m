clc;
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

global_min_x = min(x);
global_max_x = max(x);
global_min_y = min(y);
global_max_y = max(y);
global_min_z = min(z);
global_max_z = max(z);

initialWindowSize = 1000; 
stepSize = 100;           
maxLength = length(x);   

alpha_value = 2;

numWindows = floor((maxLength - initialWindowSize) / stepSize) + 1;

re_x = zeros(numWindows, 1);
re_y = zeros(numWindows, 1);
re_z = zeros(numWindows, 1);
timePoints = zeros(numWindows, 1);

for i = 1:numWindows
    currentWindowSize = initialWindowSize + (i-1)*stepSize;
    
    xWindow = x(1:currentWindowSize);
    yWindow = y(1:currentWindowSize);
    zWindow = z(1:currentWindowSize);

    re_x(i) = renyi_entropy_fixed_bins(xWindow, alpha_value, global_min_x, global_max_x);
    re_y(i) = renyi_entropy_fixed_bins(yWindow, alpha_value, global_min_y, global_max_y);
    re_z(i) = renyi_entropy_fixed_bins(zWindow, alpha_value, global_min_z, global_max_z);
    
    
    timePoints(i) = t(currentWindowSize);
end


figure;

subplot(3, 1, 1);
plot(timePoints, re_x, 'r', 'LineWidth', 1.5);
ylabel(['Rényi Entropy (α = ' num2str(alpha_value) ')']);
title('Rényi Entropy of x vs Time');
xlim([t(initialWindowSize), t(end)]);
grid on;

subplot(3, 1, 2);
plot(timePoints, re_y, 'g', 'LineWidth', 1.5);
ylabel(['Rényi Entropy (α = ' num2str(alpha_value) ')']);title('Rényi Entropy of y vs Time');
xlim([t(initialWindowSize), t(end)]);
grid on;

subplot(3, 1, 3);
plot(timePoints, re_z, 'b', 'LineWidth', 1.5);
xlabel('Time');
ylabel(['Rényi Entropy (α = ' num2str(alpha_value) ')']);title('Rényi Entropy of z vs Time');
xlim([t(initialWindowSize), t(end)]);
grid on;

set(gcf, 'Position', [100, 100, 800, 600]);

function dXdt = lorentz(~, X, a, b, c)
    x = X(1);
    y = X(2);
    z = X(3);
    
    dx = a*(y - x);
    dy = c*x - x*z - y;
    dz = x*y - b*z;
    
    dXdt = [dx; dy; dz];
end

function re = renyi_entropy_fixed_bins(data, alpha, global_min, global_max)
    
    num_bins = 100;
    bin_edges = linspace(global_min, global_max, num_bins + 1);
    
    counts = histcounts(data, bin_edges);
    probabilities = counts / sum(counts);
    
    probabilities = probabilities(probabilities > 0);
    
    if alpha == 1
        re = -sum(probabilities .* log2(probabilities));
        return;
    end
    
    re = (1 / (1 - alpha)) * log2(sum(probabilities .^ alpha));
end