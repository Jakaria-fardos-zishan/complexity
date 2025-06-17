% Parameters
a = 10; b = 8/3; c = 28;
tspan = [0 100];
dt = 0.01;
t = tspan(1):dt:tspan(2);

% Initial conditions
X0 = [1; 1; 1];

% Solve ODE
[tSol, X] = ode45(@(t, X) lorentz(t, X, a, b, c), t, X0);

% Get time series
x = X(:,1);
y = X(:,2);
N = length(x);

% Discretize signals into bins
numBins = 10;
x_disc = discretize(x, numBins);
y_disc = discretize(y, numBins);

% Prepare for TE computation
% (x_t, y_t) --> y_{t+1}
x_t   = x_disc(1:end-1);
y_t   = y_disc(1:end-1);
y_tp1 = y_disc(2:end);

% Joint histogram: p(y_{t+1}, y_t, x_t)
joint_xyz = accumarray([y_tp1 y_t x_t], 1, [numBins numBins numBins]);

% Normalize to get joint probabilities
P_xyz = joint_xyz / sum(joint_xyz(:));

% Marginal: p(y_t, x_t)
P_yx = sum(P_xyz, 1);
P_yx = squeeze(P_yx);

% Marginal: p(y_t)
P_y = sum(sum(P_xyz, 3), 1);
P_y = squeeze(P_y);

% Conditional probs
TE = 0;
for i = 1:numBins
    for j = 1:numBins
        for k = 1:numBins
            p_joint = P_xyz(i,j,k);
            if p_joint == 0
                continue;
            end
            p_yx = P_yx(j,k);
            p_y = P_y(j);
            p_y_cond = sum(P_xyz(i,j,:)); % sum over x_t for p(y_{t+1}, y_t)
            
            p1 = p_joint / p_yx; % p(y_{t+1} | y_t, x_t)
            p2 = p_y_cond / p_y; % p(y_{t+1} | y_t)

            if p1 > 0 && p2 > 0
                TE = TE + p_joint * log2(p1 / p2);
            end
        end
    end
end

fprintf('Transfer Entropy (X → Y): %.4f bits\n', TE);
