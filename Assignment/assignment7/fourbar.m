clear all
clf
clc
global lambda time i
lambda = [];time =[];
time_int = [0 12.25];

%initial condition

theta1 = 
phi1=
phi2=



theta1dot = 0;
phi2dot = 0;
phi1dot = 0;
q_vector = [theta1 phi2 phi1 theta1dot phi2dot phi1dot];


[t,q] = ode45(@four_bar,time_int,q_vector);

subplot(2,1,1)
plot(t,q(:,1)*180/pi,'linewidth',2);
hold on
plot(t,q(:,2)*180/pi,'linewidth',2);
plot(t,q(:,3)*180/pi,'linewidth',2);
legend('\theta_1','\phi_2','\phi_1'); 

for i = 1:length(t)
lambda(:,i) = lamb(q(i,:));
end
time = t;
subplot(2,1,2)
plot(time(3:length(time)),lambda(1,3:length(time)),'linewidth',2);
hold on
plot(time(3:length(time)),lambda(2,3:length(time)),'linewidth',2);
legend('\lambda_1','\lambda_2'); 



function dq = four_bar(t,q)
global lambda time i
L1 = ; L2 =; L3 =; L4 = ;
r1 = ; r2 = ; r3 = ;
I1 = ; I2 =; I3 = ;
m1 = ; m2 = ; m3 = ;
g = ;
T0 = ; k = ;
theta1 = q(1); 
phi2 = q(2);
phi1 = q(3);
theta1dot = q(4);
phi2dot = q(5);
phi1dot = q(6);

% Mass matrix
M = [];
 
 
% C matrix 
C = [];
 
% G matrix 
G = [];
 
 
 
psi = [];
    
    
    
psidot = [];
      
% Torque acting    
T = [T0-k*theta1; 0; 0];


qdot = q(4:6);
% lambda lagrange multiplier
i = i+1;
lambda(:,i) = -(inv(psi*inv(M)*psi'))*(psidot*qdot+psi*inv(M)*(T-C*qdot-G));
time(i) = t;

A = -C-psi'*inv(psi*inv(M)*psi')*(-psi*inv(M)*C+psidot);
B = T-G-psi'*inv(psi*inv(M)*psi')*psi*inv(M)*(T-G);

dq = zeros(6,1);    % a column vector 
dq(1:3) = q(4:6);
dq(4:6) = inv(M)*A*q(4:6)+inv(M)*B;
 
end


function lambda = lamb(q)

L1 = ; L2 = ; L3 = ; L4 = ;
r1 = ; r2 = L2/2; r3 = ;
I1 = ; I2 =; I3 = ;
m1 = ; m2 =; m3 = ;
g = ;
T0 = ; k =;
theta1 = q(1); 
phi2 = q(2);
phi1 = q(3);
theta1dot = q(4);
phi2dot = q(5);
phi1dot = q(6);


M = [];
 
 
 
C = [];
 
 
G = [];
 
 
 
psi = [];
    
    
    
psidot = [];
      
    
T = [T0-k*theta1; 0; 0];

% qdot = zeros(3,1);
qdot = q(4:6)';
lambda = -(inv(psi*inv(M)*psi'))*(psidot*qdot+psi*inv(M)*(T-C*qdot-G));
end
