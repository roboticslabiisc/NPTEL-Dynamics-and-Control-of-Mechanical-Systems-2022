clc
clear

a = 2; % Scalar

A1 = [1,2,3]; % Row vector

A2 = [4;5;6]; % Column vector

A3 = [1 2 8;7 5 6]; % 2x3 Matrix

A4 = [A3;A1]; % Concatenation

A5 = A4'; % Transpose

I1 = eye(4); % 4x4 identity matrix

N1 = norm(A1); % L2 norm of A1

C1 = cross(A1,A2); % Cross product

D1 = dot(A1,A2); % Dot product

B1 = det(A4); % Determinant

B2 = inv(A4); % Inverse
B5 = A4^-1;

Ra = rank(A3); % Rank of a matrix

R1 = rotx(30); % Rotation matrix for rotation about X-axis by 30 degrees
R2 = roty(45); % Rotation matrix for rotation about Y-axis by 45 degrees
R3 = rotz(60); % Rotation matrix for rotation about Z-axis by 60 degrees

R = R1*R2*R3
R4 = a*R; % Multiplying a scaler with a matrix
V1 = R*A2; % Multiplying a vector with a matrix

eul = [30 45 60]*pi/180;
rotmXYZ = eul2rotm(eul,'XYZ'); % Body-fixed euler angles to rotation matrix
eulXYZ = rotm2eul(rotmXYZ,'XYZ')*180/pi; % Rotaion matrix to body-fixed euler angles

axang = [0 1 0 pi/6]; % Axis-angle(rad)
rotmaxang = axang2rotm(axang); % Axis-angle to rotation matrix
ax_ang = rotm2axang(rotmaxang); % Rotaion matrix to axis-angle

% Similar conversions to quaternions is also present

[V,D] = eig(R2); % Eigenvalues(D) and eigenvectors(V) of R2

J = R(1,2);
K = R(2,:);
norm(R);

T = linspace(0,10,11);
E = T(2:2:end);
L = length(T);
S = size(T);
