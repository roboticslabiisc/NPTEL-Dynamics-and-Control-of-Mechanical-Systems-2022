% Solving linear equations of the form A*X=b
% Equation 1:  x - 2y = -2
% Equation 2: 7x - 3y = 19

clc
clear

A = [1,-2;7,-3];
b = [-2;19];

X = linsolve(A,b)
Y = A\b