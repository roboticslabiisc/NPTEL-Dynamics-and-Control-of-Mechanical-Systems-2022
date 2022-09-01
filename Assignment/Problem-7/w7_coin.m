clc
clear
clf

m = 0.5;%mass of the disk;
r = 0.1;%radius of the disk;
g = 9.8;%gravitational acceleration;

timestep = 0.01;
timespan = [0:timestep:10];

% Set initial conditions satisfying constraints;
% control theta1, theta2, theta3, x, y, dtheta1, dtheta2, dtheta3;

theta10 = -60*pi/180;
theta20 = 5*pi/180;
theta30 = 0*pi/180;
Xc0 = 0;
Yc0 = 0;
Zc0 = r*cos(theta20);
dtheta10 = 500*pi/180;
dtheta20 = 20*pi/180;
dtheta30 = 500*pi/180;
dXc0 = r*cos(theta10)*sin(theta20)*dtheta10 + r*sin(theta10)*cos(theta20)*dtheta20...
	+ r*cos(theta10)*dtheta30;
dYc0 = r*sin(theta10)*sin(theta20)*dtheta10 - r*cos(theta10)*cos(theta20)*dtheta20...
	+ r*sin(theta10)*dtheta30;
dZc0 = -r*sin(theta20)*dtheta20;
 
% Vector of initial conditions for ode solver;
 
options = odeset('RelTol',1e-9,'Stats','on','OutputFcn',@odeplot);
y0 = [Xc0 Yc0 Zc0 theta10 theta20 theta30 dXc0 dYc0 dZc0 dtheta10 dtheta20 dtheta30];
[t,y] = ode45(@(t,y) Lagrangian_Eqns(t,y),timespan,y0,options);

q = y(:,1:6);
dq = y(:,7:12);
Xc = y(:,1);
Yc = y(:,2);
Zc = y(:,3);
theta1 = y(:,4);
theta2 = y(:,5);
theta3 = y(:,6);

% Animate the Rolling 
% Animate = Animate_Disk(r, theta1, theta2, theta3, Xc, Yc, Zc, t,y);

figure(2)

subplot(3,2,1)
plot(t,Xc)
xlabel('Time (s)')
ylabel('X_c (m)')
grid on

subplot(3,2,3)
plot(t,Yc)
xlabel('Time (s)')
ylabel('Y_c (m)')
grid on

subplot(3,2,5)
plot(t,Zc)
% plot(t,r*cos(theta2))
xlabel('Time (s)')
ylabel('Z_c (m)')
grid on

subplot(3,2,2)
plot(t,theta1)
xlabel('Time (s)')
ylabel('\theta_1 (rad)')
grid on

subplot(3,2,4)
plot(t,theta2)
xlabel('Time (s)')
ylabel('\theta_2 (rad)')
grid on

subplot(3,2,6)
plot(t,theta3)
xlabel('Time (s)')
ylabel('\theta_3 (rad)')
grid on

KE = [];
PE = [];
for i = 1:1:length(t)
    y1 = y(i,1); y2 = y(i,2); y3 = y(i,3); y4 = y(i,4); y5 = y(i,5); y6 = y(i,6);
    y7 = y(i,7); y8 = y(i,8); y9 = y(i,9); y10 = y(i,10); y11 = y(i,11); y12 = y(i,12);
    ke = (1/2)*y7^2*m+(1/2)*y8^2*m+(1/2)*y9^2*m+(-(1/8)*y10*r^2*(cos(y5)^2-2)*m+(1/4)*y12*r^2*sin(y5)*m)*y10+(1/8)*y11^2*m*r^2+((1/4)*y10*r^2*sin(y5)*m+(1/4)*y12*m*r^2)*y12;
    PE = [PE,m*g*r*cos(y(i,5))];
    KE = [KE,ke];
    TE = KE + PE;
end

figure(3)
plot(t,KE,'r')
hold on
plot(t,PE,'b')
plot(t,TE,'k')
xlabel('Time (s)')
ylabel('Energy (J)')
% axis([0 t(end) 0 round(max(TE))])
legend({'Kinetic Energy','Potential Energy','Total Energy'},'Location','northoutside','NumColumns',3)
hold off

% calculate KE & PE;
% plot all the 3 constraint errors;

function [ ydot ] = Lagrangian_Eqns( t,y )
%LAGRANGIAN_EyNS Summary of this function goes here
%Detailed explanation goes here

%definitions
Xc = y(1);
yc = y(2);
Zc = y(3);
theta1 = y(4);
theta2 = y(5);
theta3 = y(6);
dXc = y(7);
dyc = y(8);
dZc = y(9);
dtheta1 = y(10);
dtheta2 = y(11);
dtheta3 = y(12);
q = y(1:6);
dq = y(7:12);
m = 0.5;
r = 0.1;
g = 9.8;


%generalized forces;
Q = [0; 0; 0; 0; 0; 0];

%generalized mass matrix;
M = [
	[m	0	0	0	0	0];
	[0	m	0	0	0	0];
	[0	0	m	0	0	0];
	[0	0	0	1/8*m*(-2*r^2*cos(theta2)^2+4*r^2)	0	1/2*sin(theta2)*m*r^2];
	[0	0	0	0	1/4*m*r^2	0];
	[0	0	0	1/2*sin(theta2)*m*r^2	0	1/2*m*r^2]
    ];

%generalized damping matrix;
C = [
	[0	0	0	0	0	0];
	[0	0	0	0	0	0];
	[0	0	0	0	0	0];
	[0	0	0	1/4*m*r^2*cos(theta2)*sin(theta2)*dtheta2	1/4*m*r^2*cos(theta2)*sin(theta2)*dtheta1+1/4*cos(theta2)*m*r^2*dtheta3	1/4*cos(theta2)*m*r^2*dtheta2];
	[0	0	0	-1/4*m*r^2*cos(theta2)*sin(theta2)*dtheta1-1/4*cos(theta2)*m*r^2*dtheta3	0	-1/4*cos(theta2)*m*r^2*dtheta1];
	[0	0	0	1/4*cos(theta2)*m*r^2*dtheta2	1/4*cos(theta2)*m*r^2*dtheta1	0]
	];

%generalized gravity matrix;
G = [0;0;0;0;-m*g*r*sin(theta2);0];

%constraint matrix;
psi = [
	[1	0	0	-cos(theta1)*sin(theta2)*r	-r*sin(theta1)*cos(theta2)	-r*cos(theta1)];
	[0	1	0	-sin(theta1)*sin(theta2)*r	r*cos(theta1)*cos(theta2)	-r*sin(theta1)];
	[0	0	1	0	sin(theta2)*r	0]
	];
    
%derivative of constarint matrix;   
dpsi = [
	[0	0	0	r*sin(theta1)*sin(theta2)*dtheta1-r*cos(theta2)*cos(theta1)*dtheta2	-r*cos(theta1)*dtheta1*cos(theta2)+r*sin(theta1)*sin(theta2)*dtheta2	r*sin(theta1)*dtheta1];
	[0	0	0	-r*cos(theta1)*sin(theta2)*dtheta1-r*cos(theta2)*sin(theta1)*dtheta2	-r*sin(theta1)*dtheta1*cos(theta2)-r*cos(theta1)*sin(theta2)*dtheta2	-r*cos(theta1)*dtheta1];
	[0	0	0	0	cos(theta2)*dtheta2*r	0]
	];
    
    
f = Q-C*dq-G;
lambda = -(psi/M*psi')\(dpsi*dq+psi/M*f);%Lagrange multiplier;
ddq = M\(f+psi'*lambda);%q doubledot;
ydot = [dq; ddq];

end

function [Animate] = Animate_Disk(r, theta1, theta2, theta3, Xc, Yc, Zc, t,y)
%ANIMATE Summary of this function goes here
%   Detailed explanation goes here

% r = 0.1;%radius 
phi = linspace(0, 2*pi, 50)';
c1 = r*cos(phi);
c2 = -r*sin(phi);

%Body_fixed axes {Xb,Yb,Zb}.Point A is on disk's periphery;
%Zblf is loosely fixed on disk;  

timestep = 0.01;
timespan = [0:timestep:20];

for j = 1:length(t)
  
    Rz = [cos(theta1(j)) -sin(theta1(j)) 0;sin(theta1(j)) cos(theta1(j)) 0;0 0 1];   
    Rx = [1 0 0;0 cos(theta2(j)) -sin(theta2(j));0 sin(theta2(j)) cos(theta2(j))];
    Ry = [cos(theta3(j)) 0 sin(theta3(j));0 1 0;-sin(theta3(j)) 0 cos(theta3(j))]; 
    R = Rz*Rx*Ry;
    Zblf(:,j) = (Rz*Rx)*[0;0;1];              
    Xb(:,j) = R*[1;0;0];                      
    Yb(:,j) = R*[0;1;0];                              
    Zb(:,j) = R*[0;0;1];                      
    x_on_circle(:,j) = Xc(j) + c1*Xb(1,j) + c2*Zb(1,j);     
    y_on_circle(:,j) = Yc(j) + c1*Xb(2,j) + c2*Zb(2,j);     
    z_on_circle(:,j) = Zc(j) + c1*Xb(3,j) + c2*Zb(3,j);                                                
    xP(j,1) = Xc(j) - r*Zblf(1,j);                          
    yP(j,1) = Yc(j) - r*Zblf(2,j);                          
    zP(j,1) = 0;                                            
    xA(j,1) = Xc(j) + r*Xb(1,j);                            
    yA(j,1) = Yc(j) + r*Xb(2,j);                            
    zA(j,1) = Zc(j) + r*Xb(3,j);                         
end

%Setting up the figure;
myVideo = VideoWriter('coin'); %open video file
myVideo.FrameRate = 50;  %can adjust this, 5 - 10 works well for me
open(myVideo)

figure(1);
set(gcf, 'color', 'w')
plot3(Xc(1), Yc(1), Zc(1));
xlabel('Xaxis')
set(gca, 'xdir', 'reverse')
ylabel('Yaxis')
set(gca, 'ydir', 'reverse')
zlabel('Zaxis')
axis equal
xmin = min(y(:,1)); xmax = max(y(:,1));
ymin = min(y(:,2)); ymax = max(y(:,2));
zmin = min(y(:,3)); zmax = max(y(:,3)); 
%figure;
grid on;

%view(3);
%Highlight point A on the disk's periphery;

path_of_P = line('xdata', xP(1:1), 'ydata', yP(1:1), 'zdata', zP(1:1), 'color', 'c', 'linewidth', 1);
Disk = patch('xdata', x_on_circle(:,1), 'ydata', y_on_circle(:,1), 'zdata', z_on_circle(:,1), 'facecolor', [1, 0, 1], 'linewidth', 2);
path_of_A = line('xdata', xA(1), 'ydata', yA(1), 'zdata', zA(1), 'marker', '*', 'color', 'g', 'markerfacecolor', 'g', 'linewidth', 3);

for j = 1:length(xP)
    sgtitle("Time : " + t(j) + " seconds")
    set(path_of_P, 'xdata', xP(1:j), 'ydata', yP(1:j), 'zdata', zP(1:j));
    set(Disk, 'xdata', x_on_circle(:,j), 'ydata', y_on_circle(:,j), 'zdata', z_on_circle(:,j));
    set(path_of_A, 'xdata', xA(j), 'ydata', yA(j), 'zdata', zA(j));
    drawnow 
    frame = getframe(gcf); %get frame
    writeVideo(myVideo, frame);
end
% close(myVideo)
Animate = [];
end

