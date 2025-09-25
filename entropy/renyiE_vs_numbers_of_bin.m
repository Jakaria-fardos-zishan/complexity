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

alpha_value = 2;

bin_values = 10:10:100000;

re_x = zeros(size(bin_values));
re_y = zeros(size(bin_values));
re_z = zeros(size(bin_values));


for i = 1:length(bin_values)
    num_bins = bin_values(i);
    re_x(i) = renyi_entropy_with_bins(x, alpha_value, num_bins);
    re_y(i) = renyi_entropy_with_bins(y, alpha_value, num_bins);
    re_z(i) = renyi_entropy_with_bins(z, alpha_value, num_bins);
end

figure;
plot(bin_values, re_x, 'r', 'LineWidth', 2);
hold on;
plot(bin_values, re_y, 'g', 'LineWidth', 2);
plot(bin_values, re_z, 'b', 'LineWidth', 2);
hold off;

xlabel('Number of Bins');
ylabel(['Rényi Entropy (α = ' num2str(alpha_value) ')']);
title('Rényi Entropy vs Number of Bins for Lorenz System');
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

function re = renyi_entropy_with_bins(data, alpha, num_bins)
   
    [counts, ~] = histcounts(data, num_bins);
    
    probabilities = counts / sum(counts);
    
    probabilities = probabilities(probabilities > 0);
    
    if alpha == 1
        re = -sum(probabilities .* log2(probabilities));
        return;
    end
    
    
    re = (1 / (1 - alpha)) * log2(sum(probabilities .^ alpha));
end