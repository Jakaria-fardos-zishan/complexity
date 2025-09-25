clc;
a=20;
b=8/3;
c=28;
x0=1;
y0=5;
z0=9;
time=[0,5000];

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

figure(1);
plot(x, y, 'b');
xlabel('x');
ylabel('y');
title('x vs y');
grid on;

figure(2);
plot(y, z, 'r');
xlabel('y');
ylabel('z');
title('y vs z');
grid on;

figure(3);
plot(z, x, 'g');
xlabel('z');
ylabel('x');
title('z vs x');
grid on;

figure(4);
plot(v_x, v_y, 'm');
xlabel('v_x');
ylabel('v_y');
title('v_x vs v_y');
grid on;

figure(5);
plot(v_y, v_z, 'c');
xlabel('v_y');
ylabel('v_z');
title('v_y vs v_z');
grid on;

figure(6);
plot(v_z, v_x, 'k');
xlabel('v_z');
ylabel('v_x');
title('v_z vs v_x');
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
