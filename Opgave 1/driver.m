%% Problem 1.3: EqualityQPSolver and RandomQP
H = [6,2,1; 2,5,2;1,2,4];
g = [-8;-3;-3];
A = [1,0;0,1;1,1];
b = [3;0];

[x,lambda] = EqualityQPSolver(H,g,A,b);
eig(H)
%% Problem 1.5: Random QP generator
merr=zeros(1000,2);
%Run 1000 times
for i=1:1000
    n=randi(6);
    m=randi(6);
    [H,g,A,b, xT, lambdaT] = RandomQP(n,m);
    [x,lambda] = EqualityQPSolver(H,g,A,b);

    %Calculate error
    errx = abs(xT-x);
    errlambda = abs(lambdaT-lambda);
    merr(i,1) = mean(errx);
    merr(i,2) = mean(errlambda);
end
mean(merr)
%% Problem 1.7 Sensitivities
p = zeros(5,1);
[dx, dlambda, x_approx, lambda_approx] = SensitivitiesEQP(H,g,A,b,p);

% Try with [p1, p2, p3, p4, p5] = [0, 0, 1, 0, 1] and compare with new
% corresponding b and g.
p = [0, 0, 0, 1, 1]';
gp = g+p(1:3);
bp = b+p(4:5);

% First order Taylor approximations of solution is exact
[dx, dlambda, x_approx, lambda_approx] = SensitivitiesEQP(H,g,A,b,p);
[x,lambda] = EqualityQPSolver(H,gp,A,bp);





