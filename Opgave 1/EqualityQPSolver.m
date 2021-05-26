function [x,lambda] = EqualityQPSolver(H,g,A,b)
% Solves convex quadratic program with equality constrains

%Dimensions
[s1,s2] = size(H);
[a1,a2] = size(A);

%KKT matrix
F = [H , -A; -A' , zeros(a2,a2)];


%RHS
h = [-g;-b];
eps=0;
% Use factorization, use cholesky if positive definite
if(all(eig(F) > eps))
    [L,p] = chol(F,'lower');
    x = L'\(L\h);
else % Use LDL else
    [L,D,p] = ldl(F,'vector');
    z = L'\(D\(L\h(p)));
    x = z(p);
end

%Get solution
lambda = x((s1+1):end);
x = x(1:s1);
