clc
clear
clf

t = linspace(0,2*pi,1001);

% X1 = cos(t);
% Y1 = sin(t);
% 
% plot(X1,Y1)
% grid on
% axis equal
% title('Trignometric functions')
% xlabel('X-axis')
% ylabel('Y-axis')
% axis([-1.5 2*pi -1.3 1.4])
% 
% X2 = t;
% Y2 = sin(t);
% 
% hold on
% 
% plot(X2,Y2)
% 
% legend('x^2+y^2=1','sin(x)')

X3 = sin(t);
Y3 = cos(t);
Z3 = t;

clf
plot3(X3,Y3,Z3)
axis equal
hold on

scatter3(X3(1:50:end),Y3(1:50:end),Z3(1:50:end),'rx')

clf
fimplicit(@(x,y) (x^2 + 2*y^2 -1))
axis equal