clc
clear
clf

sys1 = tf([2 5 1],[1 2 3 4]); % Defining the transfer function form
bodeplot(sys1) % Bode plot of the system

A = [0,1,0;0,0,1;-6,-11,-6];
B = [0;0;6];
C = [1 0 0];
D = 0;
sys2 = ss(A,B,C,D);  % Defining the state-space form
bodeplot(sys2) % Bode plot of the system

grid on

margin(sys2)