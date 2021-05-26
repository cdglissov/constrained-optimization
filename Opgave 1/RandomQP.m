function [H,g,A,b,x,lambda] = RandomQP(n,m)
% n is the dimensions of H
% m is the number of constraints
H = rand(n,n);
% Generate H based on random matrix
% Assure it's symmetric
H = 0.5*(H+H');
% Assure it's positive definite
H = H+n*eye(n);
% Generate A
A = rand(n,m);
%Generate true solution
x = rand(n,1);
% Construct b
b = A'*x;
% Define Lambda
lambda = rand(m,1);
% Construct g
g = A*lambda-H*x;