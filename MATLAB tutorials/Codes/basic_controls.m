clc
clear
clf

% Defining state-space matrices
A1 = [0 1 0;0 0 1;-1 -5 -6];
B1 = [0;0;1];
C1 = [0 0 1];
D1 = 0;
[num1,den1] = ss2tf(A1,B1,C1,D1); % Convert state-space to transfer function
sys1 = ss(A1,B1,C1,D1);


% Defining transfer function form
num2 = [0 1 0 0];
den2 = [1 6 5 1];
[A2,B2,C2,D2] = tf2ss(num2,den2); % Convert transfer function to state-space form
sys2 = tf(num2,den2);

[P,V] = jordan(A1); % Jordan fom

syms t
eAt = expm(A2*t);

[z,p,k] = zpkdata(sys2); % Obtaining the zeros(z), poles(p) and gain(k) of the system


t = 0:0.01:25; % Defining the time for obtaining the response
ic = [1;0;0]; % Defining the initial conditions of state vector
initial(sys1,ic,t) % Evolution of output(y)


step(sys2) % Step response
impulse(sys1) % Impulse response

syms a t
f = exp(-a*t); % Defining a function of time
L = laplace(f); % Obtaining laplace transform of a function (f)
I = ilaplace(L); % Obtaining the inverse laplace