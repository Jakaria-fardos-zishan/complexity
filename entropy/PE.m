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

initialWindowSize = 1000; 
stepSize = 100;          
maxLength = length(x);  


m = 3;
tau = 1;


numWindows = floor((maxLength - initialWindowSize) / stepSize) + 1;

pe_x = zeros(numWindows, 1);
pe_y = zeros(numWindows, 1);
pe_z = zeros(numWindows, 1);
timePoints = zeros(numWindows, 1);

for i = 1:numWindows
   
    currentWindowSize = initialWindowSize + (i-1)*stepSize;
    
    xWindow = x(1:currentWindowSize);
    yWindow = y(1:currentWindowSize);
    zWindow = z(1:currentWindowSize);
    
  
    pe_x(i) = permutation_entropy(xWindow, m, tau);
    pe_y(i) = permutation_entropy(yWindow, m, tau);
    pe_z(i) = permutation_entropy(zWindow, m, tau);
    
    timePoints(i) = t(currentWindowSize);
end

figure;

subplot(3, 1, 1);
plot(timePoints, pe_x, 'r', 'LineWidth', 1.5);
ylabel('Permutation Entropy');
title('Permutation Entropy of x vs Time');
xlim([t(initialWindowSize), t(end)]);
grid on;

subplot(3, 1, 2);
plot(timePoints, pe_y, 'g', 'LineWidth', 1.5);
ylabel('Permutation Entropy');
title('Permutation Entropy of y vs Time');
xlim([t(initialWindowSize), t(end)]);
grid on;

subplot(3, 1, 3);
plot(timePoints, pe_z, 'b', 'LineWidth', 1.5);
xlabel('Time');
ylabel('Permutation Entropy');
title('Permutation Entropy of z vs Time');
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

function pe = permutation_entropy(data, m, tau)
    n = length(data);
    permutations = perms(1:m);
    n_patterns = size(permutations, 1);
    count = zeros(1, n_patterns);

    for i = 1:(n - (m - 1) * tau)
        
        delay_vector = data(i:tau:i + (m - 1) * tau);
        
      
        [~, idx] = sort(delay_vector);
        
       
        for j = 1:n_patterns
            if isequal(idx', permutations(j, :))
                count(j) = count(j) + 1;
                break;
            end
        end
    end
    
 
    total = sum(count);
    p = count(count > 0) / total;
   
    pe = -sum(p .* log2(p));
end