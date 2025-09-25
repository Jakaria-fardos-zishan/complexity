clc;
a = 10;
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

bin_range = 10:10:100000; 
se_bins_x = zeros(size(bin_range));
se_bins_y = zeros(size(bin_range));
se_bins_z = zeros(size(bin_range));

for i = 1:length(bin_range)
    num_bins = bin_range(i);
    se_bins_x(i) = shannon_entropy_with_bins(x, num_bins);
    se_bins_y(i) = shannon_entropy_with_bins(y, num_bins);
    se_bins_z(i) = shannon_entropy_with_bins(z, num_bins);
end

figure(4);
plot(bin_range, se_bins_x, 'r', 'LineWidth', 2);
hold on;
plot(bin_range, se_bins_y, 'g', 'LineWidth', 2);
plot(bin_range, se_bins_z, 'b', 'LineWidth', 2);
hold off;
xlabel('Number of Bins');
ylabel('Shannon Entropy');
title('Shannon Entropy vs Number of Bins');
legend('x', 'y', 'z');
grid on;

function dXdt = lorentz(~, X, a, b, c)
    x = X(1);
    y = X(2);
    z = X(3);
    
    dx = a*(y - x);
    dy = c*x - x*z - y;
    dz = x*y - b*z;
    
    dXdt = [dx; dy; dz];
end

function se = shannon_entropy_with_bins(data, num_bins)
    [counts, ~] = histcounts(data, num_bins);
    probabilities = counts / sum(counts);
    probabilities = probabilities(probabilities > 0);
    se = -sum(probabilities .* log2(probabilities));
end