function [x, lambda] = EQQP(H, g, A, b)

[n,m] = size(A);

% Define KKT System
KKT = [H -A; -A' zeros(m,m)];
rhs = - [g; b];

% LDL Factorize
[L, D, p] = ldl(KKT, 'vector');

% Back substitute
xlambda(p,1) = L'\(D\(L\rhs(p,1)));

% Extract x and lambda
x = xlambda(1:n,1);
lambda = xlambda(n+1:n+m,1);
