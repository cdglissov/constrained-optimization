function [K,h] = KKTSystem(n,u,d,sparsity)

% Construct H,g,A,b
[H, g, A, b] = ConstructEqQP(n,u,d,sparsity);
%Get dimensions of A
[~,a2] = size(A);
% Construct KKT matrix
K=[H , -A; -A' , zeros(a2,a2)];
%Create RHS
h = [-g;-b];


