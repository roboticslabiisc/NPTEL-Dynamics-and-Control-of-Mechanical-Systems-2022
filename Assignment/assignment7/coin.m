clc
clear
clf

% parameters given in the question
m = %mass of the disk;
r = %radius of the disk;
g = %gravitational acceleration;

timestep = 0.01;
timespan = [0:timestep:10];

% Set initial conditions satisfying constraints;
% control theta1, theta2, theta3, x, y, dtheta1, dtheta2, dtheta3;


%=========% initial conditions%==================%

theta10 =       % θ1(0)
theta20 =       % 
theta30 =       % 
Xc0 =           % x(0)
Yc0 = 
Zc0 = 
dtheta10 =      % θ˙1(0)  theta_dot(0)
dtheta20 = 
dtheta30 = 
dXc0 =          % x˙(0)  x_dot(0)
dYc0 =   
dZc0 = 
 


% Vector of initial conditions for ode solver;
 
options = odeset('RelTol',1e-9,'Stats','on','OutputFcn',@odeplot);
y0 = [Xc0 Yc0 Zc0 theta10 theta20 theta30 dXc0 dYc0 dZc0 dtheta10 dtheta20 dtheta30]; % initial condition
[t,y] = ode45(@(t,y) Lagrangian_Eqns(t,y),timespan,y0,options); % solver 

q = y(:,1:6);
dq = y(:,7:12);
Xc = y(:,1);
Yc = y(:,2);
Zc = y(:,3);
theta1 = y(:,4);
theta2 = y(:,5);
theta3 = y(:,6);



figure(2)
% This is one plot example 
% you can make subplots please checkout sub plots

plot(t,Xc)
xlabel('Time (s)')
ylabel('X_c (m)')
grid on




%Kinetic and potential energy

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
%Plots for energies
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
%LAGRANGIAN_EQNS Summary of this function goes here
%Detailed explanation goes here

%definitions please fill the respective values from y few values are already
%filled to give tou hint.
Xc = y(1);
yc = 
Zc = 
theta1 = y(4);
theta2 =
theta3 = 
dXc = y(7);
dyc = 
dZc = 
dtheta1 = y(10);
dtheta2 =
dtheta3 = 
q = y(1:6);
dq = 
m = 
r = 
g = 


%generalized forces;
Q = [0; 0; 0; 0; 0; 0];

%generalized mass matrix;
M = []; 

%generalized damping matrix;
C = [];  

%generalized gravity matrix;
G = []; 

%constraint matrix;
psi = [];
    
%derivative of constarint matrix;   
dpsi = [];
    
    
f = Q-C*dq-G;
lambda = -(psi/M*psi')\(dpsi*dq+psi/M*f);%Lagrange multiplier;
ddq = M\(f+psi'*lambda);%q doubledot;
ydot = [dq; ddq];

end

%plots
