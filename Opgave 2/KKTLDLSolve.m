function [x,lambda] = KKTLDLSolve(n,H,g,A,b,K,h)
z = zeros(2*n+1,1);
% Use LDL factorization
[L,D,p] = ldl(K,'lower','vector');
z(p) = L'\(D\(L\h(p)));
%Get solution and Langrangian multipliers
x=z(1:(n+1));
lambda=z((n+2):end);