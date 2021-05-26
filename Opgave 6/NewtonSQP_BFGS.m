function [x, stat] = NewtonSQP_BFGS(fundfun,consfun,x0,y0,varargin)
% y is Lagrange multiplier

maxit = 100*length(x0);
tol   = 1e-6;


stat.converged = false; % converged
stat.nfun = 0; % number of function calls
stat.iter = 0; % number of iterations

% Initial iteration
x = x0;
it = 0;
B = eye(numel(x0));
[f,df,d2f] = feval(fundfun,x,varargin{:});
[c,dc,d2c] = feval(consfun,x,varargin{:});
[~,y0] = EQQP(B,df,dc,c);
y=y0;
%Compute dL(x,y_new)
dL_2 = df - dc*y;
converged = (norm(dL_2,'inf') <= tol && norm(c) <= tol);
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
stat.d2C = d2c;
stat.Errc = norm(c, "inf");
stat.ErrL = norm(y, "inf");     
    
    
    while ~converged && (it < maxit)
        % updating the iteration number
        it = it + 1;
        
        
        % Solve equality constriant
        [p,y] = EQQP(B,df,dc,-c);
        
        
        % Take step
        x = x + p;
        
        
        % Function evaluation
        [f,df,d2f] = feval(fundfun,x,varargin{:});
        [c,dc,d2c] = feval(consfun,x,varargin{:});
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
        
        %Weighting
        r = theta*q + (1-theta)*B*p;
        
        dL_2=dL;
        %Approximate hessian
        B = B + (r*r')/(p'*r) - ((B*p)*(B*p)')/(p'*B*p);
        
        converged = (norm(dL,'inf') <= tol && norm(c) <= tol);
        stat.X   = [stat.X x];
        stat.Y   = [stat.Y y];
        stat.F   = [stat.F f];
        stat.C   = [stat.C c];
        stat.B   = [stat.B B];
        stat.dF  = [stat.dF df];
        stat.dC  = [stat.dC dc];
        stat.d2F = [stat.d2F d2f];
        stat.d2C = [stat.d2C d2c];
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