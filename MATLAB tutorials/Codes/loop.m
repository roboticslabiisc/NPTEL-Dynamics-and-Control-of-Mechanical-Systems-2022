% Collect all the numbers in 0-50 which are not divisible by 5 into groups:
% 1. Which are also divisible by 3 or 7.
% 2. The remaining numbers

clc
clear

t = 0:1:50;

X = [];
Y = [];
for i = 1:1:length(t)
    c = t(i);
    if mod(c,5) ~= 0
        if mod(c,3) == 0
            X = [X,c];
        elseif mod(c,7) == 0
            X = [X,c];
        else
            Y = [Y,c];
        end
    end
end

X
Y