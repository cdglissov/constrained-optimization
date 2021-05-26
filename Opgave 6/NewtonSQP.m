function [x, stat] = NewtonSQP(fundfun,consfun,x0,y0,varargin)
% y is Lagrange multiplier

maxit = 100*length(x0);
tol   = 1e-5;


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
dL_2 = df - dc*y;
converged = (norm(dL_2,'inf') <= tol && norm(c) <= tol);
stat.nfun = 2;

% Store data for plotting
stat.X = x;
stat.Y = y;
stat.F = f;
stat.C = c;
stat.dF = df;
stat.dC = dc;
stat.d2F = d2f;
stat.d2C = d2c;
stat.Errc = norm(c, "inf");
stat.ErrL = norm(dL_2, "inf");  
    
    
    while ~converged && (it < maxit)
        % updating the iteration number
        it = it + 1;
        
        % Computing the Hessian
        H = d2f;
        for i = 1:length(y0)
            H = H - y(i)*d2c(:,:,i);
        end
        
        % Solve equality constriant
        [p,y] = EQQP(H,df,dc,-c);
        
        % Take step
        x = x + p;
        
        % Function evaluation
        [f,df,d2f] = feval(fundfun,x,varargin{:});
        [c,dc,d2c] = feval(consfun,x,varargin{:});
        stat.nfun = stat.nfun + 2;
        
        dL = df - dc*y;
        
        converged = (norm(dL,'inf') <= tol && norm(c,'inf') <= tol);
        stat.X   = [stat.X x];
        stat.Y   = [stat.Y y];
        stat.F   = [stat.F f];
        stat.C   = [stat.C c];
        stat.dF  = [stat.dF df];
        stat.dC  = [stat.dC dc];
        stat.d2F = [stat.d2F d2f];
        stat.d2C = [stat.d2C d2c];
        stat.Errc = [stat.Errc norm(c, "inf")];
        stat.ErrL = [stat.ErrL norm(dL,"inf")];
    end
    

if ~converged
    x = [];
end
stat.converged = converged;
stat.iter = it;