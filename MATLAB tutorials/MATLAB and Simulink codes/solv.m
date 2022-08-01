clc
clear
clf

t = 0:0.1:4;

f1 = @(x)(x.^3 - 6*x.^2 + 11.*x - 6);

y = f1(t);
plot(t,y)
grid on

x0 = 2.1;

X = fsolve(f1,x0)