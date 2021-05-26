function [x, stat] = SQP_TRUST_REG(ObjFun1,ConFun1,x0,y0)
%SQP TRUST REGION
debug=0;
maxit = 100*length(x0);
setQP = optimoptions('quadprog','display','off');
setLP = optimoptions('linprog','display','off');
stat.converged = false; % converged
stat.nfun = 0; % number of function calls
stat.iter = 0; % number of iterations

%Initialise
x = x0;
[f,df,d2f] = feval(ObjFun1,x);
[c,dc] = feval(ConFun1,x);
stat.nfun = 2;
it = 0;
%Set hyper params
mu1=1e-5; %1e-4 also seems to work great
tol   = 1e-6;
dkmax=2;
dk = 1;  %set this for different solution, best for 1?
eps1 = 0.4;
eps2 = 0.4;
eta = 0.2; %0.1-0.3
gamma = 0.4; %0.3?
converged = 0;
B = eye(numel(x0));
m = size(c,1);
n = size(x,1);
t=ones(n,1);

% Store data for plotting
stat.X = x;
stat.F = f;
stat.C = c;
stat.B = B;
stat.dF = df;
stat.dC = dc;
stat.d2F = d2f;     
stat.Errc = norm(c, "inf");
stat.ErrL = norm(y0, "inf"); 
%mk is a piecewise linear model (18.56)
mk = @(c,dc,p) sum(max(0,-(c + dc'*p)));

while ~converged && (it < maxit)
    % updating the iteration number
    it = it + 1;
    %Used to solve 18.50
    B2 = [B, zeros(n,m); zeros(m,n), zeros(m,m)];
    %Add slack to g
    g = [df; mu1 * t];
    %t adds extra dimensions, multiply with - to swap inequality
    A = -[dc', eye(n)];
    %inequality limits
    b=c;
    if(debug==1)
        A
        B2
        g
    end
    % ||p||_inf < delta_k -> can go from -dk to dk
    % slack variables can go from t>=0 -> 0 to inf
    lb = [-dk*t; 0*t];
    ub = [dk*t; Inf*t];

    % Solve 18.50, to get slack t and p
    [pslack,~,~,~,l] = quadprog(B2,g,A,b,[],[],lb,ub,[],setQP);
    %get solution for p
    l_hat=l.ineqlin(1:2);
    p = pslack(1:2);
    converged =(norm(p,'inf') < tol);
    if(converged==1)
        stat.X   = [stat.X x];
        stat.F   = [stat.F f];
        stat.C   = [stat.C c];
        stat.B   = [stat.B B];
        stat.dF  = [stat.dF df];
        stat.dC  = [stat.dC dc];
        stat.d2F = [stat.d2F d2f];
        break;
    end
    m_k=mk(c,dc,p);
    %%% PENALTY UPDATE AND STEP COMPUTATION %%%
    if abs(m_k) <= tol
        mu = mu1;
    else
        p_lp = linprog([zeros(n,1); ones(m,1)], A, b, [], [],lb,ub,[], ...
                    setLP);
        pinf = p_lp(1:2);
        mkinf=mk(c,dc,pinf);
        if abs(mkinf) <= tol 
            %new constrained inequality problem
            mkpos=m_k;
            mu = mu1;
            %keep doing until below is satisfied
            while(abs(mkpos)>tol)
                %multiply mu by a factor, thats how we get muk>muk-1
                mu = mu*5;
                %solve the quadprog to make p(mu^pos) and mk(mu^+) satisfy
                [muslack,~,~,~,~] = quadprog(B2,[df; mu * t],A,b,[],[], ...
                                lb,ub,[],setQP);
                pmupos = muslack(1:2);
                mkpos = mk(c,dc,pmupos);
                if(debug==1)
                    pmupos
                    mkpos
                    mu
                end
            end
        else
            mk0 = mk(c,dc,zeros(2,1));
            mkpos = m_k;
            mu = mu1;
            %Keep doing until below is satisfied
            while((mk0-mkpos) < (eps1*(mk0-mkinf)))
                %scale
                mu = mu*5;
                [muslack,~,~,~,~] = quadprog(B2,[df; mu * t],A,b,[],[], ...
                                lb,ub,[],setQP);
                 %solve the quadprog to make p(mu^pos) and mk(mu^+) satisfy
                 pmupos = muslack(1:2);
                 mkpos = mk(c,dc,pmupos);
                if(debug==1)
                    pmupos
                    mkpos
                    mu
                end
            end
        end
    end
    %Quadprog again...
    [muslack,~,~,~,~] = quadprog(B2,[df; mu * t],A,b,[],[], ...
                                lb,ub,[],setQP);
    pmupos = muslack(1:2);
    %find these variables piece wise linear
    mk0 = mk(c,dc,zeros(2,1));
    mkpos = mk(c,dc,pmupos);
    %Calculate (18.57)
    q0 = f + mu*mk0;
    qp = f + df'*pmupos + 0.5*pmupos'*B*pmupos + mu*mkpos;
    if(debug==1)
        q0
        qp
    end
    %need to satisfy this scale mu until satisfied
    while((q0-qp) < eps2*mu*(mk0-mkpos))
        mu=mu*5;
        [muslack,~,~,~,~] = quadprog(B2,[df; mu * t],A,b,[],[], ...
                        lb,ub,[],setQP);
        %update
        pmupos=muslack(1:2);
        mkpos = mk(c,dc,pmupos);
        q0 = f + mu*mk0;
        qp = f + df'*pmupos + 0.5*pmupos'*B*pmupos + mu*mkpos;
    end
    
    %set mu_k = mupos and p=p(mu+)
    mu1=mu;
    p=pmupos;
    %%%% COMPLETE %%%%
    [f1,~,~] = feval(ObjFun1,x+p);
    [c1,~] = feval(ConFun1,x+p);
    stat.nfun = stat.nfun + 2;
    %%% MERIT FUNCTIONS %%%%
    % (18.48) to judge our step
    phi1 = f + mu1*sum(max(-c, 0));
    phi1p = f1 + mu1*sum(max(-c1, 0));
    if(debug==1)
        phi1
        phi1p
    end
    rhok = (phi1 - phi1p) / (q0 - qp);
    if rhok > eta
        %Update parameters
        dLold = df - dc*l_hat;
        x = x + p;

        % Function evaluation
        [f,df,d2f] = feval(ObjFun1,x);
        [c,dc] = feval(ConFun1,x);
        stat.nfun = stat.nfun + 2;

        dLnew = df - dc*l_hat;
        %compute q
        q = dLnew - dLold;

        %Update Hessian by modified damped BFGS
        if ( p'*q >= 0.2*p'*B*p)
            theta = 1;
        else
            theta = ( 0.8*p'*B*p ) /(p'*B*p - p'*q );
        end

        r = theta*q + (1-theta)*B*p;
        %approximate hessian
        B = B + (r*r')/(p'*r) - ((B*p)*(B*p)')/(p'*B*p);
        
        %update radius of circle
        dk = min(dk * 1.5,dkmax);
    else
        %update radius of circle
        dk = min(gamma*norm(p, 'inf'), dkmax);
    end
    if(debug==1)
        dk
    end
    %Save stuff
    converged =(norm(p,'inf') < tol);
    stat.X   = [stat.X x];
    stat.F   = [stat.F f];
    stat.C   = [stat.C c];
    stat.B   = [stat.B B];
    stat.dF  = [stat.dF df];
    stat.dC  = [stat.dC dc];
    stat.d2F = [stat.d2F d2f];
    stat.Errc = norm(c, "inf");
    stat.ErrL = norm(dLnew, "inf");
end
    
stat.converged = converged;
stat.iter = it;
if ~converged
    x = [];
end
stat.converged = converged;
stat.iter = it;