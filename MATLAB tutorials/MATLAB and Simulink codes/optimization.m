clc
clear

x0 = [0;0]; % Initial conditions

[sol,fval] = fmincon(@optifun,x0,[],[],[],[],[],[],@const)

function F = optifun(x)
    F = -((x(1))^2 + (x(2))^2); % Defining the optimization function
end

function [c,ceq] = const(x)
    ceq = [(x(1))^2 + (x(2))^2 - 1]; % Defining the equality constraints
    c = []; % Defining the in-equality constraints
end