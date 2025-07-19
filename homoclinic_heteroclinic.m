% Lorenz System Parameter Space Scanner (a,b,c ∈ [0,30])
% Homoclinic = Red, Heteroclinic = Blue

% Set parameter range and resolution
param_range = [0, 30];
step = 2;  % Reduced step for better resolution
tol = 1e-2;  % Tolerance for orbit detection

% Prepare parameter grid
[a_vals, b_vals, c_vals] = meshgrid(param_range(1):step:param_range(2), ...
                                    param_range(1):step:param_range(2), ...
                                    param_range(1):step:param_range(2));
homoclinic_points = [];
heteroclinic_points = [];

% Counters
total_points = numel(a_vals);
fprintf('Scanning %d parameter sets...\n', total_points);

% Main scanning loop
for i = 1:total_points
    a = a_vals(i);
    b = b_vals(i);
    c = c_vals(i);
    
    % Progress indicator
    if mod(i, 1000) == 0
        fprintf('Progress: %.1f%%\n', 100*i/total_points);
    end
    
    % Find equilibrium points
    try
        syms x y z
        eq1 = a*(y - x) == 0;
        eq2 = c*x - x*z - y == 0;
        eq3 = x*y - b*z == 0;
        eq_points = solve([eq1, eq2, eq3], [x, y, z]);
        equilibria = double([eq_points.x, eq_points.y, eq_points.z]);
    catch
        continue;  % Skip if solution fails
    end
    
    % Identify saddles (points with mixed stability)
    saddles = [];
    for j = 1:size(equilibria, 1)
        pt = equilibria(j,:);
        if any(isnan(pt)) || any(isinf(pt))
            continue
        end
        
        % Jacobian matrix
        J = [-a, a, 0; 
             c - pt(3), -1, -pt(1); 
             pt(2), pt(1), -b];
        
        % Compute eigenvalues
        eigs = eig(J);
        if isreal(eigs) && any(eigs > 0) && any(eigs < 0)
            saddles = [saddles; pt];
        end
    end
    
    % Skip if not enough saddles
    if size(saddles, 1) < 1
        continue;
    end
    
    % Test trajectory from first saddle
    saddle = saddles(1,:);
    J_s = [-a, a, 0; c - saddle(3), -1, -saddle(1); saddle(2), saddle(1), -b];
    [V, D] = eig(J_s);
    eigenvalues = diag(D);
    unstable_idx = find(real(eigenvalues) > 0, 1);
    
    if isempty(unstable_idx)
        continue;
    end
    
    % Perturb along unstable direction
    perturb = real(V(:, unstable_idx)) * 1e-5;
    X0 = saddle' + perturb;
    
    % Integrate trajectory
    try
        [~, X] = ode45(@(t,X)lorentz(t,X,a,b,c), [0 50], X0, odeset('RelTol',1e-6));
    catch
        continue;
    end
    
    % Check trajectory endpoints
    start_point = X(1,:);
    end_point = X(end,:);
    
    % Find closest saddle to endpoints
    [~, start_idx] = min(vecnorm(saddles - start_point, 2, 2));
    [~, end_idx] = min(vecnorm(saddles - end_point, 2, 2));
    
    start_dist = norm(start_point - saddles(start_idx,:));
    end_dist = norm(end_point - saddles(end_idx,:));
    
    % Classify if within tolerance
    if start_dist < tol && end_dist < tol
        if start_idx == end_idx
            homoclinic_points = [homoclinic_points; a, b, c];
        else
            heteroclinic_points = [heteroclinic_points; a, b, c];
        end
    end
end

% Visualization
figure('Position', [100, 100, 1200, 900]);
hold on;
grid on;

% Plot homoclinic parameters
if ~isempty(homoclinic_points)
    scatter3(homoclinic_points(:,1), homoclinic_points(:,2), homoclinic_points(:,3), ...
             40, 'r', 'filled', 'MarkerFaceAlpha', 0.7, 'DisplayName', 'Homoclinic');
end

% Plot heteroclinic parameters
if ~isempty(heteroclinic_points)
    scatter3(heteroclinic_points(:,1), heteroclinic_points(:,2), heteroclinic_points(:,3), ...
             40, 'b', 'filled', 'MarkerFaceAlpha', 0.7, 'DisplayName', 'Heteroclinic');
end

% Mark classic Lorenz parameters
plot3(10, 8/3, 28, 'gp', 'MarkerSize', 15, 'MarkerFaceColor', 'g', 'DisplayName', 'Classic Chaotic');
plot3(10, 8/3, 13.926, 'mo', 'MarkerSize', 10, 'DisplayName', 'Homoclinic (ρ≈13.93)');
plot3(10, 8/3, 24.74, 'co', 'MarkerSize', 10, 'DisplayName', 'Heteroclinic (ρ≈24.74)');

% Labels and title
xlabel('Parameter a (σ)');
ylabel('Parameter b (β)');
zlabel('Parameter c (ρ)');
title('Lorenz System: Special Orbits in Parameter Space (a,b,c ∈ [0,30])');
legend('Location', 'best');
view(135, 30);  % Optimal viewing angle
axis equal;
box on;
set(gca, 'FontSize', 12);

% Save memory by clearing large arrays
clear a_vals b_vals c_vals

fprintf('Scan complete. Found %d homoclinic and %d heteroclinic parameter sets.\n', ...
        size(homoclinic_points,1), size(heteroclinic_points,1));

% Lorenz system equations
function dX = lorentz(~, X, a, b, c)
    x = X(1); y = X(2); z = X(3);
    dX = [a*(y - x); 
          c*x - x*z - y; 
          x*y - b*z];
end