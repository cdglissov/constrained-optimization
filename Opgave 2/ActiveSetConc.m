function [x, W, Lambda] = ActiveSetConc(H, g, A, b)

N_c = length(b);

% Guess working set
W = dec2bin(2^N_c-1:-1:0)-'0'
P = randperm(1:2^N_c)

for i = P
    W = find(W(i,:));
    if isempty(W)
        %algo
    else
        %algo
    end
    
end
