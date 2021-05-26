function [x,lambda] = KKTNSSolve(n,H,g,A,b,K,h)
%Perform null space factorization
[Q,Rbar] = qr(A);
m1 = size(Rbar,2);
Q1 = Q(:,1:m1); Q2 = Q(:,m1+1:n+1); R = Rbar(1:m1,1:m1);
xy = R'\b;
xz = (Q2'*H*Q2)\(-Q2'*(H*Q1*xy+g));
x = Q1*xy+Q2*xz;
lamb = R\(Q1'*(H*x+g));
z = vertcat(x,lamb);
%Get solution and Langrangian multipliers
x=z(1:(n+1));
lambda=z((n+2):end);