function [x, stat] = SQP_BFGS_ineq(ObjFun1,ConFun1,x0,y0)
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
B = eye(numel(x0));

[f,df,d2f] = feval(ObjFun1,x);
[c,dc] = feval(ConFun1,x);
y=y0;
%Compute dL(x,y_new)
dL_2 = df - dc*y;
converged = 0;
stat.nfun = 2;

% Store data for plotting
stat.X = x;
stat.Y = y;
stat.F = f;
stat.C = c;
stat.B = B;
stat.dF = df;
stat.dC = dc;
stat.d2F = d2f;
stat.Errc = norm(c, "inf");
stat.ErrL = norm(y, "inf");     
    
    
    while ~converged && (it < maxit)
        % updating the iteration number
        it = it + 1;
        
        % Solve equality constriant
        [p,~,~,~,l] = quadprog(B,df,-dc',c,[],[],[],[],[],setQP);
        y=l.ineqlin;
        
        % Take step
        x = x + p;
        
        % Function evaluation
        [f,df,d2f] = feval(ObjFun1,x);
        [c,dc] = feval(ConFun1,x);
        stat.nfun = stat.nfun + 2;
        
        dL = df - dc*y;

        %compute q
        q = dL - dL_2;
        
        %Update Hessian by modified BFGS
        if ( p'*q >= 0.2*p'*B*p)
            theta = 1;
        else
            theta = ( 0.8*p'*B*p ) /(p'*B*p - p'*q );
        end
        
        r = theta*q + (1-theta)*B*p;
        dL_2=dL;
        %Approximate Hessian
        B = B + (r*r')/(p'*r) - ((B*p)*(B*p)')/(p'*B*p);
        
        
        
        converged =(norm(p,'inf') < tol);
        stat.X   = [stat.X x];
        stat.Y   = [stat.Y y];
        stat.F   = [stat.F f];
        stat.C   = [stat.C c];
        stat.B   = [stat.B B];
        stat.dF  = [stat.dF df];
        stat.dC  = [stat.dC dc];
        stat.d2F = [stat.d2F d2f];
        stat.Errc = [stat.Errc norm(c, "inf")];
        stat.ErrL = [stat.ErrL norm(dL,"inf")]; 
    end
    
stat.converged = converged;
stat.iter = it;
if ~converged
    x = [];
end
stat.converged = converged;
stat.iter = it;