function [f, df, d2f] = ObjFun(x)

f1 = x(1)^2 + x(2) - 11;
f2 = x(1) + x(2)^2 -7;

% Objective function
f = f1.^2 + f2.^2;

% Gradient of f
df = zeros(2,1);
df(1,1) = 4*x(1)*f1 + 2*f2;
df(2,1) = 2*f1 + 4*x(2)*f2;

% Hessian of f
d2f = zeros(2,2);
d2f(1,1) = 4*f1 + 8*x(1)^2 + 2;
d2f(1,2) = 4*(x(1)+x(2));
d2f(2,1) = d2f(1,2);
d2f(2,2) = 4*f2 + 8*x(2)^2 +2;