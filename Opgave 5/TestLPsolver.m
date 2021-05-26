%
% Test LP solver
%

n = 200;
m = 50;

state = 1000;
%rand('state',state);

A = randn(m,n);

x = zeros(n,1);
x(1:m,1) = abs(rand(m,1));

lambda = zeros(n,1);
lambda(m+1:n,1) = abs(rand(n-m,1));

mu = rand(m,1);

g = A'*mu + lambda;
b = A*x;

%g, A, b

[xlp,info,mulp,lambdalp,iter] = LPippd(g,A,b,ones(n,1));

iter=iter

if info
    X = max(abs(xlp-x))
    M = max(abs(mulp-mu))
    L = max(abs(lambdalp-lambda))
end
