function [dx, dlambda, x_approx, lambda_approx] = SensitivitiesEQP(H,g,A,b,p)
% Computing solution for p = 0, in order to approximate x(p) and lambda(p)
[x0,lambda0] = EqualityQPSolver(H,g,A,b);

 % s1 number of var, s2 number of constrains
[s1,s2] = size(A);

% Determined parameters in Q1.6
dxc = A;
Wxx = H;
z = zeros(s2,s2);

%Setup sensitivity matrix
K = [Wxx -dxc; -dxc' z];

%Sensitivity is calculated
Kinv = -inv(K);

% Sensitivity of x(p) and sensitivity of lambda(p)
dx = Kinv(:,1:s1);
dlambda = Kinv(:,s1+(1:s2));

% 1st order Taylor approximation to the solutions
x_approx = x0 + dx'*p;
lambda_approx = lambda0 + dlambda'*p;





