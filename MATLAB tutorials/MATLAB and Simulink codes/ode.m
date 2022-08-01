clc
clear
clf

global m c k f

m = 1;
c = 1;
k = 1;
f = 0;

time = 0:0.1:10;
y0 = [1;0];

[t,y] = ode45(@eqn,time,y0);

plot(t,y(:,1))
hold on
plot(t,y(:,2))
xlabel('Time(s)')
legend('x_1','x_2')

function F = eqn(t,x)
    global m c k f
    
    F1 = x(2);
    F2 = (f - c*x(2) - k*x(1))/m;
    F = [F1;F2];
end