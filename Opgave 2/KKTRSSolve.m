function [x,lambda] = KKTRSSolve(n,H,g,A,b,K,h)
%Perform range-space factorization using a cholesky factorization
L = chol(H,'lower');
v = L'\(L\g);
Ha = A'*inv(L)'*inv(L)*A;
La = chol(Ha, 'lower');
lamb = La'\(La\(b+A'*v));
x = L'\(L\(A*lamb-g));
z = vertcat(x,lamb);
%Get solution and Langrangian multipliers
x=z(1:(n+1));
lambda=z((n+2):end);