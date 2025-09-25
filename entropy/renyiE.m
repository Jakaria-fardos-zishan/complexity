% Clear workspace and command window
clear;
clc;

% Lorenz system parameters
a = 10;
b = 8/3;
c = 28;
initial_conditions = [1; 5; 9];
time_span = [0, 5000];

% Solve Lorenz system
[t, XYZ] = ode45(@(t, X) lorentz(t, X, a, b, c), time_span, initial_conditions);
x = XYZ(:, 1);
y = XYZ(:, 2);
z = XYZ(:, 3);

% Sliding window parameters
windowSize = 1000;
stepSize = 100;
numWindows = floor((length(x) - windowSize) / stepSize) + 1;
alpha_value = 2; % Rényi entropy parameter

% Precompute global bin edges for each variable using all data
numBins = 100;
xEdges = linspace(min(x), max(x), numBins + 1);
yEdges = linspace(min(y), max(y), numBins + 1);
zEdges = linspace(min(z), max(z), numBins + 1);

% Preallocate arrays for Rényi entropy and time midpoints
re_x = zeros(numWindows, 1);
re_y = zeros(numWindows, 1);
re_z = zeros(numWindows, 1);
timeMidpoints = zeros(numWindows, 1);

% Calculate Rényi entropy for each window
for i = 1:numWindows
    startIdx = (i-1)*stepSize + 1;
    endIdx = startIdx + windowSize - 1;
    timeWindow = t(startIdx:endIdx);
    
    % Compute entropy with fixed bin edges
    re_x(i) = renyi_entropy_fixed_bins(x(startIdx:endIdx), alpha_value, xEdges);
    re_y(i) = renyi_entropy_fixed_bins(y(startIdx:endIdx), alpha_value, yEdges);
    re_z(i) = renyi_entropy_fixed_bins(z(startIdx:endIdx), alpha_value, zEdges);
    
    timeMidpoints(i) = mean(timeWindow);
end

% Plot for x
figure;
plot(timeMidpoints, re_x, 'r', 'LineWidth', 1.5);
xlabel('Time');
ylabel(['Rényi Entropy (\alpha = ', num2str(alpha_value), ')']);
title('Rényi Entropy of x vs Time');
xlim([50, 5000]);
grid on;

% Plot for y
figure;
plot(timeMidpoints, re_y, 'g', 'LineWidth', 1.5);
xlabel('Time');
ylabel(['Rényi Entropy (\alpha = ', num2str(alpha_value), ')']);
title('Rényi Entropy of y vs Time');
xlim([50, 5000]);
grid on;

% Plot for z
figure;
plot(timeMidpoints, re_z, 'b', 'LineWidth', 1.5);
xlabel('Time');
ylabel(['Rényi Entropy (\alpha = ', num2str(alpha_value), ')']);
title('Rényi Entropy of z vs Time');
xlim([50, 5000]);
grid on;

% Lorenz system differential equations
function dXdt = lorentz(~, X, a, b, c)
    x = X(1);
    y = X(2);
    z = X(3);
    dXdt = [a*(y - x); 
            c*x - x*z - y; 
            x*y - b*z];
end

% Rényi entropy with fixed bin edges
function re = renyi_entropy_fixed_bins(data, alpha, edges)
    counts = histcounts(data, edges);
    probabilities = counts(counts > 0) / sum(counts);
    if alpha == 1
        re = -sum(probabilities .* log2(probabilities));
    else
        re = (1 / (1 - alpha)) * log2(sum(probabilities .^ alpha));
    end
end