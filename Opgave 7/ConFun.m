function [c,dc] = ConFun(x)

% Constraints
c1 = (x(1)+2)^2 - x(2);
c2 = -4*x(1)+10*x(2);

% Gradient of c1
dc1 = zeros(2,1);
dc1(1,1) = 2*x(1)+4;
dc1(2,1) = -1;

% Gradient of c2
dc2 = zeros(2,1);
dc2(1,1) = -4;
dc2(2,1) = 10;

% Combine constraints
c = [c1; c2];
dc = [dc1, dc2];