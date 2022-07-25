clc
clear

z = [-1 -5];
p = [-4 -2 0 3];
k = 1;

sys = zpk(z,p,k);

[num,den] = zpkdata(sys,'v');

[A,B,C,D] = tf2ss(num,den);
K = place(A,B,[-4.5 -1.5 -1])
eig(A)

% [num,den] = ss2tf(A,-B.*K',C,D)


eig(A-B*K)
