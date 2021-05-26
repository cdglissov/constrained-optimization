function [c,dc,d2c]=nlpcon(x)


n=5;
m=3;

% Constraint functions
c = zeros(m,1);
c(1,1) = sum(x.*x) - 10;
c(2,1) = x(2)*x(3) - 5*x(4)*x(5);
c(3,1) = x(1)^3 + x(2)^3 + 1;


% Gradient functions
dc = zeros(n,m);

dc(:,1) = 2*x;

dc(2,2) = x(3);
dc(3,2) = x(2);
dc(4,2) = -5*x(5);
dc(5,2) = -5*x(4);

dc(1,3) = 3*x(1)^2;
dc(2,3) = 3*x(2)^2;

% Hessians
d2c = zeros(n,n,m);

d2c(:,:,1) = 2*eye(n,n);

d2c(3,2,2) = 1.0;
d2c(2,3,2) = 1.0;
d2c(5,4,2) = -5.0;
d2c(4,5,2) = -5.0;

d2c(1,1,3) = 6*x(1);
d2c(2,2,3) = 6*x(2);
