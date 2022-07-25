clc
clear

syms phi a b c
k = [1,0,0];
R = rotmat(k,phi);

R_xyz = rotmat([1,0,0],a)*rotmat([0,1,0],b)*rotmat([0,0,1],c)


function R = rotmat(k,phi)
    K = skew(k);
    R = eye(3) + K*sin(phi) + (1-cos(phi))*K^2;
end

function W = skew(k)
    W = [0,-k(3),k(2);
         k(3),0,-k(1);
         -k(2),k(1),0];
end