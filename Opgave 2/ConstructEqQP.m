function [H, g, A, b] = ConstructEqQP(n,ub,d0,sparsity)
% Construct g
g = -ones(n+1,1)*ub;
e1 = ones(n,1);
e2 = zeros(n,1);
% Construct H and A sparse or dense
if sparsity == 1
    H = sparse(eye(n+1,n+1));
    A = spdiags([e1 -e1 e2],[-1 0 1],n,n+1); %sparse
elseif sparsity == 0
    H = eye(n+1,n+1);
    A = spdiags([e1 -e1 e2],[-1 0 1],n,n+1);
    A = full(A); % dense
else
    disp('Choose 0 or 1 sparsity');
end
% Finalizing A
A(1,n) = 1;
A(n,n+1) = -1;
A = A';
% Constructing b
b = zeros(n,1);
b(1) = -d0;