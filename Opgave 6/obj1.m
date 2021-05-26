function [g,dg,d2g]=obj1(x)


n = size(x,1);

% Function
t = prod(x);
g = exp(t);

% Gradient
dt = zeros(n,1);
for i=1:n
    k = setdiff(1:n,i);
    dt(i,1) = prod(x(k,1));
end
dg = dt.*g;

% Hessian
d2g = (dg)*dt';
for j=1:n
    for i=1:j-1
        k = setdiff(1:n,[i j]);
        d2g(i,j) = d2g(i,j) + prod(x(k,1))*g;
    end
    
    for i=j+1:n
        k = setdiff(1:n,[i j]);
        d2g(i,j) = d2g(i,j) + prod(x(k,1))*g;
    end
end
