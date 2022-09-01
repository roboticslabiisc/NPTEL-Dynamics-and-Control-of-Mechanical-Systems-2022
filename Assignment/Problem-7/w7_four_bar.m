clc
clear
clf

global m1 m2 m3 I1 I2 I3 L1 L2 L3 L0 r1 r2 r3 tau0 K theta10 phi10 g

g = 9.81;

% m1 = 20;
% m2 = 10;
% m3 = 10;
% I1 = 9;
% I2 = 0.1;
% I3 = 0.1;
% 
% L1 = 1.241;
% L2 = 1.2;
% L3 = 1.2;
% L0 = 1.241;
% r1 = 1.2;
% r2 = 0.6;
% r3 = 0.6;
% 
% tau0 = 2.0;
% K = 0.1;

% Book
m1 = 20.15;
m2 = 8.25;
m3 = 8.25;
I1 = 9.6;
I2 = 0.06;
I3 = 0.06;

L1 = 1.241;
L2 = 1.2;
L3 = 1.2;
L0 = 1.241;
r1 = 1.2;
r2 = 0.6;
r3 = 0.6;

tau0 = 1.96;
K = 0.1;

tspan = 0:0.01:10;
theta10 = 0.01;
phi10 = 0.0102;
phi20 = 6.2698; %6.2698

y0 = [theta10;phi20;phi10;0;0;0];

[t,y] = ode45(@eom,tspan,y0);
q1 = y(:,1); q2 = y(:,2); q3 = y(:,3); q4 = y(:,4); q5 = y(:,5); q6 = y(:,6);

L = [];
for i = 1:1:length(t)
    L = [L,lambda([q1,q2,q3,q4,q5,q6])];
end

figure(1)
plot(t,q1,t,q2,t,q3)
legend('\theta_1','\phi_2','\phi_1')

figure(2)
plot(t,L(1,:),t,L(2,:))

function Q = eom(t,q)
    global m1 m2 m3 I1 I2 I3 L1 L2 L3 r1 r2 r3 tau0 K g
    q1 = q(1); q2 = q(2); q3 = q(3); q4 = q(4); q5 = q(5); q6 = q(6);
    dq = [q4;q5;q6];
    
    M = [0.2e1 * cos(q2) * L1 * r2 * m2 + (L1 ^ 2 + r2 ^ 2) * m2 + r1 ^ 2 * m1 + I1 + I2 cos(q2) * L1 * r2 * m2 + m2 * r2 ^ 2 + I2 0; cos(q2) * L1 * r2 * m2 + m2 * r2 ^ 2 + I2 m2 * r2 ^ 2 + I2 0; 0 0 m3 * r3 ^ 2 + I3;];
    C = [-sin(q2) * L1 * r2 * m2 * q5 -sin(q2) * L1 * m2 * r2 * (q4 + q5) 0; sin(q2) * L1 * r2 * m2 * q4 0 0; 0 0 0;];
    G = [m1 * g * r1 * cos(q1) + m2 * g * (L1 * cos(q1) + r2 * cos(q1 + q2)) cos(q1 + q2) * g * m2 * r2 m3 * g * r3 * cos(q3)]';
    psi = [-L1 * sin(q1) - L2 * sin(q1 + q2) -L2 * sin(q1 + q2) L3 * sin(q3); L1 * cos(q1) + L2 * cos(q1 + q2) L2 * cos(q1 + q2) -L3 * cos(q3);];
    dpsi = [-L1 * q4 * cos(q1) - L2 * (q4 + q5) * cos(q1 + q2) -L2 * (q4 + q5) * cos(q1 + q2) L3 * q6 * cos(q3); -L1 * q4 * sin(q1) - L2 * (q4 + q5) * sin(q1 + q2) -L2 * (q4 + q5) * sin(q1 + q2) L3 * q6 * sin(q3);];
    tau = [tau0-K*q1;0;0];
    f = tau - C*dq - G;
%     lambda = -(psi/M*psi')\(dpsi*dq+psi/M*f);%Lagrange multiplier;
%     ddq = M\(f+psi'*lambda);%q doubledot;
    ddq = M^-1 *(f-(psi')*(psi*M^-1 *psi')^-1 *(psi*M^-1 *f + dpsi*dq));
    Q=[dq;ddq];
end

function Q = lambda(q)
    global m1 m2 m3 I1 I2 I3 L1 L2 L3 r1 r2 r3 tau0 K g
    q1 = q(1); q2 = q(2); q3 = q(3); q4 = q(4); q5 = q(5); q6 = q(6);
    dq = [q4;q5;q6];
    
    M = [0.2e1 * cos(q2) * L1 * r2 * m2 + (L1 ^ 2 + r2 ^ 2) * m2 + r1 ^ 2 * m1 + I1 + I2 cos(q2) * L1 * r2 * m2 + m2 * r2 ^ 2 + I2 0; cos(q2) * L1 * r2 * m2 + m2 * r2 ^ 2 + I2 m2 * r2 ^ 2 + I2 0; 0 0 m3 * r3 ^ 2 + I3;];
    C = [-sin(q2) * L1 * r2 * m2 * q5 -sin(q2) * L1 * m2 * r2 * (q4 + q5) 0; sin(q2) * L1 * r2 * m2 * q4 0 0; 0 0 0;];
    G = [m1 * g * r1 * cos(q1) + m2 * g * (L1 * cos(q1) + r2 * cos(q1 + q2)) cos(q1 + q2) * g * m2 * r2 m3 * g * r3 * cos(q3)]';
    psi = [-L1 * sin(q1) - L2 * sin(q1 + q2) -L2 * sin(q1 + q2) L3 * sin(q3); L1 * cos(q1) + L2 * cos(q1 + q2) L2 * cos(q1 + q2) -L3 * cos(q3);];
    dpsi = [-L1 * q4 * cos(q1) - L2 * (q4 + q5) * cos(q1 + q2) -L2 * (q4 + q5) * cos(q1 + q2) L3 * q6 * cos(q3); -L1 * q4 * sin(q1) - L2 * (q4 + q5) * sin(q1 + q2) -L2 * (q4 + q5) * sin(q1 + q2) L3 * q6 * sin(q3);];
    tau = [tau0-K*q1;0;0];
    f = tau - C*dq - G;
    Q = -(psi/M*psi')\(dpsi*dq+psi/M*f);%Lagrange multiplier;
end

function F = loop(x)
    global L1 L2 L3 L0 theta10 phi10
    F = (L1*cos(theta10) - L0 - L3*cos(phi10))^2 + (L1*sin(theta10) - L3*sin(phi10))^2 - L2^2;
end
