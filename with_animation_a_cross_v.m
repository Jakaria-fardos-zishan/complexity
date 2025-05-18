clc;
a=10;
b=8/3;
c=28;
x0=1;
y0=5;
z0=9;



time=[0,500];

[t,XYZ]=ode45(@(t,X) lorentz(t,X,a,b,c), time, [x0;y0;z0] );

x=XYZ(:,1);
y=XYZ(:,2);
z=XYZ(:,3);

N = length(t);

dXdt = zeros(length(t), 3);
for i = 1:length(t)
    dXdt(i,:) = lorentz(t(i), XYZ(i,:)', a, b, c)';
end 

v_x = dXdt(:,1);
v_y = dXdt(:,2);
v_z = dXdt(:,3);
v_vector = [v_x, v_y, v_z];

dt = [0; diff(t)]; 
%v_x = gradient(x)./dt;
%v_y = gradient(y)./dt;
%v_z = gradient(z)./dt;
%v_vector = [v_x, v_y, v_z];

 
a_x = a*(v_y - v_x);

a_y = c*v_x - v_x.*z - v_z.*x - v_y;
a_z = v_x.*y + v_y.*x - b*v_z;

a_vector = [a_x, a_y, a_z];

dot_va = dot(v_vector, a_vector, 2);

cross_va = cross(a_vector, v_vector, 2);

vdota_vcrossa =  cross_va .* dot_va;

% --- Set up Figure 1: 3D Trajectory ---
figure(1);
h_traj = plot3(nan, nan, nan, 'b-', 'LineWidth', 1.5); hold on;
h_point = plot3(nan, nan, nan, 'ro', 'MarkerFaceColor', 'r');
xlabel('cross\_va\_x'); ylabel('cross\_va\_y'); zlabel('cross\_va\_z');
title('Animated 3D Trajectory of cross(v, a)');
grid on; axis equal

% --- Set up Figure 2: 3 Subplots ---
figure(2);
ax1 = subplot(3,1,1); h1 = plot(ax1, t(1), vdota_vcrossa(1,1), 'b');
ylabel('vdot\_vcrossa\_x'); title('vdot\_vcrossa\_x vs Time'); grid on;
ax2 = subplot(3,1,2); h2 = plot(ax2, t(1), vdota_vcrossa(1,2), 'r');
ylabel('vdot\_vcrossa\_y'); title('vdot\_vcrossa\_y vs Time'); grid on;
ax3 = subplot(3,1,3); h3 = plot(ax3, t(1), vdota_vcrossa(1,3), 'g');
ylabel('vdot\_vcrossa\_z'); xlabel('Time'); title('vdot\_vcrossa\_z vs Time'); grid on;

xlim(ax1, [min(t) max(t)]); xlim(ax2, [min(t) max(t)]); xlim(ax3, [min(t) max(t)]);
ylim(ax1, [min(vdota_vcrossa(:,1)) max(vdota_vcrossa(:,1))]);
ylim(ax2, [min(vdota_vcrossa(:,2)) max(vdota_vcrossa(:,2))]);
ylim(ax3, [min(vdota_vcrossa(:,3)) max(vdota_vcrossa(:,3))]);

% --- Animation Loop ---
for k = 2:length(t)
    % Update 3D Trajectory
    set(h_traj, 'XData', cross_va(1:k,1), 'YData', cross_va(1:k,2), 'ZData', cross_va(1:k,3));
    set(h_point, 'XData', cross_va(k,1), 'YData', cross_va(k,2), 'ZData', cross_va(k,3));
    figure(1); drawnow limitrate;
    
    % Update Subplots
    set(h1, 'XData', t(1:k), 'YData', vdota_vcrossa(1:k,1));
    set(h2, 'XData', t(1:k), 'YData', vdota_vcrossa(1:k,2));
    set(h3, 'XData', t(1:k), 'YData', vdota_vcrossa(1:k,3));
    figure(2); drawnow limitrate;
    
    pause(0.01); % Adjust animation speed as needed
end

function dXdt= lorentz (~,X,a,b,c)
x=X(1);
y=X(2);
z=X(3);


dx= a*(y - x);
dy = c*x - x*z - y;
dz = x*y - b*z;

dXdt=[dx;dy;dz];
end


