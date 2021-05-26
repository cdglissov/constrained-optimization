function [h,dh,d2h]=obj2(x)

n = 5;
t1 = x(1)*x(1);
t2 = x(2)*x(2);
t = t1*x(1) + t2*x(2) + 1;

% Function
h = 0.5*t*t;

% Gradient
dh = zeros(n,1);
dh(1,1) = t*3*t1;
dh(2,1) = t*3*t2;

% Hessian
d2h = zeros(n,n);
d2h(1,1) = 9*t1*t1 + t*6*x(1);
d2h(2,1) = 9*t1*t2;
d2h(1,2) = 9*t1*t2;
d2h(2,2) = 9*t2*t2 + t*6*x(2);
