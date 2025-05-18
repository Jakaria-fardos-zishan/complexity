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
   
%sum=zeros(N,3);

figure;
subplot(7,1,1);
plot(t, cross_va(:,1), 'r');
ylabel('cross\_va\_x');
title('cross\_va\_x vs Time');
grid on;

subplot(7,1,2);
plot(t, cross_va(:,2), 'g');
ylabel('cross\_va\_y');
title('cross\_va\_y vs Time');
grid on;

subplot(7,1,3);
plot(t, cross_va(:,3), 'b');
ylabel('cross\_va\_z');
xlabel('Time');
title('cross\_va\_z vs Time');
grid on;

subplot(7,1,4);
plot(t, dot_va );
ylabel('cross\_va\_z');
xlabel('Time');
title('dot_va vs Time');
grid on;


subplot(7,1,5);
plot(t, vdota_vcrossa(:,1));
ylabel('vdot\_vcrossa\_x');
xlabel('Time');
title('vdot\_vcrossa\_x vs Time');
grid on;


subplot(7,1,6);
plot(t, vdota_vcrossa(:,2));
ylabel('vdot\_vcrossa\_y');
xlabel('Time');
title('vdot\_vcrossa\_y vs Time');
grid on;

subplot(7,1,7);
plot(t, vdota_vcrossa(:,3));
ylabel('vdot\_vcrossa\_z');
xlabel('Time');
title('vdot\_vcrossa\_z vs Time');
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
