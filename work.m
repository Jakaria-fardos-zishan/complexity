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

dx = [0; diff(x)];  
dy = [0; diff(y)];   
dz = [0; diff(z)]; 
 
a_x = a*(v_y - v_x);

a_y = c*v_x - v_x.*z - v_z.*x - v_y;
a_z = v_x.*y + v_y.*x - b*v_z;

a_vector = [a_x, a_y, a_z];

dW = a_x.*dx + a_y.*dy + a_z.*dz;

w = cumsum(dW);

figure;
plot(t, w, 'LineWidth', 1.5);
xlabel('Time');
ylabel('Work');
title('Work for Unit Mass');
grid on;

function dXdt= lorentz (~,X,a,b,c)
x=X(1);
y=X(2);
z=X(3);


dx= a*(y - x);
dy = c*x - x*z - y;
dz = x*y - b*z;

dXdt=[dx;dy;dz];
end

