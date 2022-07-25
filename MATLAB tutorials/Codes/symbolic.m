clc
clear

syms a b c da db dc t s

Rx = [1,0,0; 0,cos(a),-sin(a); 0,sin(a),cos(a)];
Ry = [cos(b),0,sin(b); 0,1,0; -sin(b),0,cos(b)];
Rz = [cos(c),-sin(c),0; sin(c),cos(c),0; 0,0,1];

R = Rx*Ry*Rz;

R_dot = diff(Rx,a)*Ry*Rz*da + Rx*diff(Ry,b)*Rz*db + Rx*Ry*diff(Rz,c)*dc;
Omega_R = simplify(R_dot*transpose(R));
omega_s = [Omega_R(3,2);Omega_R(1,3);Omega_R(2,1)];
omega_b = simplify(transpose(R)*omega_s);

A = [0,1,0;0,0,1;1,-3,3];
E1 = expm(A*t); % Exponential of a matrix

K = (s*eye(3) - A)^-1;
E2 = ilaplace(K); % Exponential of a matrix using inverse laplace

simplify(E1-E2);