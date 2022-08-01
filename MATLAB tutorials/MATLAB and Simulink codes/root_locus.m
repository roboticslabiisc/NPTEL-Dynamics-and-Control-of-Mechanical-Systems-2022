clc
clear
clf

sys1 = tf([2 5 1],[1 2 3 4]); % Defining the transfer function form
rlocus(sys1) % Root locus of the system
grid on
axis equal

A = [0,1,0;0,0,1;-6,-11,-6];
B = [0;0;6];
C = [1 0 0];
D = 0;
sys2 = ss(A,B,C,D); % Defining the state-space form
% rlocus(sys2)
[r,k] = rlocus(sys2); % Obtaining the location of roots(r) as the gain(k) changes

% Animating the movement of roots as k changes
for i = 1:1:length(k)
    rlocus(sys2)
    hold on
    scatter(real(r(1,i)),imag(r(1,i)),'b*')
    scatter(real(r(2,i)),imag(r(2,i)),'g*')
    scatter(real(r(3,i)),imag(r(3,i)),'r*')
    text(-4,4,strcat('K=',num2str(k(i))))
    pause(0.01)
    clf
end