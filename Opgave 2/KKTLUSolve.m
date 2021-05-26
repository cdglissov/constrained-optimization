function [x,lambda] = KKTLUSolve(n,H,g,A,b,K,h)
% Make LU Factorization
[L,U,p] = lu(K,'vector');
z(p) = U\(L\h(p));
%Get solution and Langrangian multipliers
x=z(1:(n+1));
lambda=z((n+2):end);