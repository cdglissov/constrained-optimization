function [H,g,A,b,C,d] = ConstructMarkowitz(R, rf);
%Construct minimizer
if rf==0
    H=[2.3, 0.93, 0.62, 0.74, -0.23;
        0.93, 1.4, 0.22, 0.56, 0.26;
        0.62, 0.22, 1.8, 0.78, -0.27;
        0.74, 0.56, 0.78, 3.4, -0.56;
        -0.23, 0.26, -0.27, -0.56, 2.6];
    A = [15.1, 12.5, 14.7, 9.02, 17.68;
    1 1 1 1 1]';
elseif rf==1
   H = [ 2.30 0.93 0.62 0.74 -0.23 0;
        0.93 1.40 0.22 0.56 0.26 0;
        0.62 0.22 1.80 0.78 -0.27 0;
        0.74 0.56 0.78 3.40 -0.56 0;
        -0.23 0.26 -0.27 -0.56 2.60 0;
        0 0 0 0 0 0];
    A = [15.10 12.50 14.70 9.02 17.68 2; 1 1 1 1 1 1]';
end
dim2 = size(H,1);
g=zeros(dim2,1);
%Construct equality constraint
dim1 = size(A,1);
b = [R;1];
%Construct inequality
C=eye(dim1);
d = zeros(dim1,1);
end