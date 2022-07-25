clc
clear

% Example 1
A = [0,1,0;0,0,1;-6,-11,-6];
B = [0;0;6];
C = [1 0 0];
Q = ctrb(A,B);
rank(Q)
P = obsv(A,C);
rank(P)

% Example 2
A = [0,1,0;0,0,1;0,-2,-3];
B = [0;0;1];
C = [3 4 1];
Q = ctrb(A,B);
rank(Q)
P = obsv(A,C);
rank(P)
 
% Example 3
A = [0,-0.4;1,-1.3];
B = [0.8;1];
C = [0 1];
Q = ctrb(A,B)
rank(Q)
P = obsv(A,C)
rank(P)