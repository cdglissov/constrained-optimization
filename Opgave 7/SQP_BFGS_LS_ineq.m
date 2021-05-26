function [x, stat] = SQP_BFGS_LS_ineq2(ObjFun1,ConFun1,x0,y0)
% y is Lagrange multiplier

maxit = 100*length(x0);
tol   = 1e-6;
setQP = optimoptions('quadprog','display','off');

stat.converged = false; % converged
stat.nfun = 0; % number of function calls
stat.iter = 0; % number of iterations


% Initial iteration
x = x0;
it = 0;
eta = 0.4; % (0,0.5)
tau = 0.8; % (0,1)
rho = 0.8; % (0,1)
B = eye(numel(x0));

[f,df,d2f] = feval(ObjFun1,x);
[c,dc] = feval(ConFun1,x);
l_k=y0;
%Compute dL(x,y_new)
converged =0;
stat.nfun = 2;

% Store data for plotting
stat.X = x;
stat.Y = l_k;
stat.F = f;
stat.C = c;
stat.B = B;
stat.dF = df;
stat.dC = dc;
stat.d2F = d2f;   
    

while ~converged && (it < maxit)
    % updating the iteration number
    it = it + 1;

    % Solve equality constriant
    [p,~,~,~,l] = quadprog(B,df,-dc',c,[],[],[],[],[],setQP);


    l_hat=l.ineqlin;

    %P_lambda:
    pl = l_hat - l_k;

    %%%%%  Back Tracking Line Search
    % Choose mu k
    mu = max(((df'*p+(1/2)*p'*B*p)/((1-rho)*norm(c,1))),0); %new

    %Set alpha
    alpha = 1;

    [f1,~,~] = feval(ObjFun1,x+alpha*p);
    [c1,~] = feval(ConFun1,x+alpha*p);
    %Calculate phi, l1 merit function, see week 8, slide 23:
    phi_1 = f1+mu*norm(c1,1); %perhaps update c
    %directional derivative D is found p. 542
    phi_2 = f + mu*norm(c,1) + eta*alpha*(df'*p-mu*norm(c,1));
    %Scale alpha until criterion is satisfied
    while phi_1 > phi_2
        alpha = tau*alpha;
        stat.nfun = stat.nfun + 2;
        [f1,~,~] = feval(ObjFun1,x+alpha*p);
        [c1,~] = feval(ConFun1,x+alpha*p);
        phi_1 = f1+mu*norm(c1,1);
        phi_2 = f + mu*norm(c,1) + eta*alpha*(df'*p-mu*norm(c,1));
    end

    %compute gradient lagrangian
    dL_old = df - dc*l_k;

    %update lagrangian
    l_k = l_k + alpha*pl;

    %Take a step
    x=x+alpha*p;

     % Function evaluation
    [f,df,d2f] = feval(ObjFun1,x);
    [c,dc] = feval(ConFun1,x);
    stat.nfun = stat.nfun + 2;

    dL_new = df - dc*l_k;


    %compute q
    q = dL_new - dL_old;
    p=alpha*p;
    %Update Hessian by modified BFGS
    if ( p'*q >= 0.2*p'*B*p)
        theta = 1;
    else
        theta = ( 0.8*p'*B*p ) /(p'*B*p - p'*q );
    end

    r = theta*q + (1-theta)*B*p;

    B = B + (r*r')/(p'*r) - ((B*p)*(B*p)')/(p'*B*p);


    converged =(norm(p,'inf') < tol);
    stat.X   = [stat.X x];
    stat.Y   = [stat.Y l_k];
    stat.F   = [stat.F f];
    stat.C   = [stat.C c];
    stat.B   = [stat.B B];
    stat.dF  = [stat.dF df];
    stat.dC  = [stat.dC dc];
    stat.d2F = [stat.d2F d2f];
end

stat.converged = converged;
stat.iter = it;
if ~converged
    x = [];
end
stat.converged = converged;
stat.iter = it;