function [x, stat] = NewtonSQP_lineSearch(fundfun,consfun,x0,y0,varargin)
% l_opt is Lagrange multiplier
rng(2)
maxit = 100*length(x0);
tol   = 1e-6;
k = 0;
reg1 = 0.9;
stat.converged = false; % converged
stat.nfun = 0; % number of function calls
stat.iter = 0; % number of iterations

% Initial iteration
x = x0;

B = eye(numel(x0));
[f,df,d2f] = feval(fundfun,x,varargin{:});
[c,dc,d2c] = feval(consfun,x,varargin{:});
[~,y0] = EQQP(B,df,dc,c);
y=y0;
lambda = abs(y);
%Initialise gradient lagrangian
dL_2 = df - dc*y;
converged = (norm(dL_2,'inf') <= tol && norm(c, 'inf') <= tol);
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
stat.ErrL = norm(dL_2, "inf");  

while ~converged && (k < maxit)

    % updating the iteration number
    k = k + 1;

    % Solve equality constriant
    [p,l_opt] = EQQP(B,df,dc,-c);

    pl = l_opt - y;

    %%%%%  Line Search
    stop = 0;
    alpha = 1;
    %Powell
    lambda = max(abs(l_opt),1/2*(lambda+abs(l_opt))); 
    cls = f + lambda'*abs(c);
    b = df'*p - lambda'*abs(c); 

    %START LINE SEARCH
    while stop == 0
        xk = x + alpha*p;

        [f,~,~] = feval(fundfun,xk,varargin{:});
        [c,~,~] = feval(consfun,xk,varargin{:});
        stat.nfun = stat.nfun + 2;
        phi_alpha = f + lambda'*abs(c);
        if phi_alpha <= (cls + (1-reg1)*b*alpha ) 
            stop = 1;
        else
            %Find good alpha for step length
            a = (phi_alpha - (cls + b*alpha))/(alpha^2);
            alpha_min = -b/(2*a);
            alpha = min( reg1*alpha,max(alpha_min,(1-reg1)*alpha));               
        end
    end

    x = xk; %new

    %update lagrangian
    y = y + alpha*pl;
%       % Function evaluation
    [f,df,d2f] = feval(fundfun,x,varargin{:});
    [c,dc,d2c] = feval(consfun,x,varargin{:});
%       stat.nfun = stat.nfun + 2;

    %%%% BFGS approximation
    dL = df - dc*y;

%       %compute q
    q = dL - dL_2;
%         
%       %Update Hessian by modified BFGS
    if ( p'*q > 0.2*p'*B*p)
         theta = 1;
    else
         theta = ( 0.8*p'*B*p ) /(p'*B*p - p'*q );
    end
%         
    r = theta*q + (1-theta)*B*p;
%Hessian approximation         
    B = B + (r*(r'))/(p'*r) - ((B*p)*((B*p)'))/(p'*B*p);
    dL_2 = dL; 

    converged = (norm(dL,'inf') <= tol && norm(c, 'inf') <= tol);
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
    stat.ErrL = [stat.ErrL norm(dL_2,inf)];
end

if ~converged
    x = [];
end
stat.converged = converged;
stat.iter = k;
