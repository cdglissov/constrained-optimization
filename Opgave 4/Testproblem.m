% week 7 
H = [2 0; 0 2];
g = [-2; -5];
A = [1,0;1,0]';
b = [1;1];
C = [1 -2; -1 -2; -1 2; 1 0; 0 1]';
d = [-2; -6; -2; 0; 0];
x = [1;0];
y = ones(2,1);
z = ones(5,1);
s = ones(5,1);

[x,info,mu,iter] = PDPCIP(H,g,A,C,b,d,x,y,z,s)