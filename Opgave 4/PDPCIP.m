function [x,y,z,s,info,mu,iter] = PDPCIP(H,g,A,C,b,d,x,y,z,s)

%%
[n,m] = size(A);
% n variables
% m equality contrains
nc = length(y); % Number of equality constrains
mc = length(z); % Number of inequality consrains


maxit = 100;
tolL = 1.0e-5;
tolA = 1.0e-5;
tolC = 1.0e-5;
tolSZ = 1.0e-5;
tolmu = 1.0e-5;

eta = 0.995;

% Residuals
rL = H*x+g-A*y-C*z;
rA = b-A'*x;
rC = s+d-C'*x;
rSZ = diag(s)*diag(z)*ones(length(z),1);
mu = (z'*s)/mc;

% Converged
Converged = (norm(rL,inf) <= tolL) && ...
            (norm(rA,inf) <= tolA) && ...
            (norm(rC,inf) <= tolC) && ...
            (norm(rSZ,inf) <= tolSZ) && ...
            (abs(mu) <= tolmu);

%%        
iter = 0;

while ~Converged && (iter<maxit)
    iter = iter+1;
    
    % ====================================================================
    % Form and Factorize Matrix
    % ====================================================================
    zdivs = z./s;
    sdivz = s./z;
    H1 = H + C*diag(zdivs)*C';
    K = [H1 -A; -A' zeros(m,m)];
    [L,D,p] = ldl(K,'vector');      % factorization
    
    % ====================================================================
    % Affine Step
    % ====================================================================
    % Solve
    xyaff = zeros(m+nc,1);
    temp = rSZ./z;
    temp2 = diag(zdivs)*(rC-temp);
    rL1 = rL - C*temp2;
    rhs = -[rL1; rA];
    xyaff(p) = L'\(D\(L\rhs(p)));       % back substitution
    xaff = xyaff(1:n);
    
    % Find z and s
    zaff = -diag(zdivs)*C'*xaff+temp2;
    saff = -temp-diag(sdivz)*zaff;
    
    % Step length
    izaff = find(zaff < 0.0);
    alpha1 = min([1.0; -z(izaff,1)./zaff(izaff,1)]);
    
    isaff = find(saff < 0.0);
    beta1 = min([1.0; -s(isaff,1)./saff(isaff,1)]);
    
    alpha = min(alpha1,beta1);
    
    % ====================================================================
    % Center Parameter and duality gap
    % ====================================================================
    muaff = (z+alpha*zaff)'*(s+alpha*saff)/mc;
    sigma = (muaff/mu)*(muaff/mu)*(muaff/mu);

    % ====================================================================
    % Affine-Centering-Correction Direction
    % ====================================================================
    rSZ1 = rSZ + saff.*zaff - sigma*mu*ones(mc,1);
    
    % Solve
    dxy = zeros(m+nc,1);
    temp = rSZ1./z;
    temp2 = diag(zdivs)*(rC-temp);
    rL1 = rL - C*temp2;
    rhs = -[rL1; rA];
    
    dxy(p) = L'\(D\(L\rhs(p)));       % back substitution
    dx = dxy(1:n);
    dy = dxy(n+1:end);
    
    % Find z and s
    dz = -diag(zdivs)*C'*dx+temp2;
    ds = -temp-diag(sdivz)*dz;
    
    % Step length
    iz = find(dz < 0.0);
    alpha1 = min([1.0; -z(iz,1)./dz(iz,1)]);
    
    is = find(ds < 0.0);
    beta1 = min([1.0; -s(is,1)./ds(is,1)]);
    
    alpha = min(alpha1,beta1);
    
    % ====================================================================
    % Update iteration
    % ====================================================================
    alpha1 = eta*alpha;
    x = x + alpha1*dx;
    y = y + alpha1*dy;
    z = z + alpha1*dz;
    s = s + alpha1*ds;
    
    % Residuals
    rL = H*x+g-A*y-C*z;
    rA = b-A'*x;
    rC = s+d-C'*x;
    rSZ = diag(s)*diag(z)*ones(length(z),1);
    mu = (z'*s)/mc;
    
    % Converged
    Converged = (norm(rL,inf) <= tolL) && ...
            (norm(rA,inf) <= tolA) && ...
            (norm(rC,inf) <= tolC) && ...
            (norm(rSZ,inf) <= tolSZ) && ...
            (abs(mu) <= tolmu);
end

%%

% Return solution
info = Converged;
if ~Converged
    x=[];
    mu=[];
    lambda=[];
end