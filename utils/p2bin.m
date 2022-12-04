function [T, rind, tind] = p2bin(R, b)
% documentation
% [T, rind, tind] = P2Bin(R, b)
% inputs:
% R: reference (= input) image
% b: binning factor (any power of 2)
% outputs:
% T: 
% rind: linear index to access reference
% tind: linear index to access template

[M, N] = size(R);
if or( mod(M, 2) ~= 0, mod(N, 2) ~= 0)
    error('#rows and #columns of reference need to be powers of 2!')
end

if and(mod(b, 2) ~= 0, b~=1)
    error('binning factor needs to be a power of 2')
end

if b~=1
    T = R;
    for k = 1:log2(b)
        T = bin2(T);
    end
    
    tind = 1:numel(T);
    rind = reshape(1:(M*N), [M, N]);
    rind = mat2cell(rind, b*ones(M/b,1), b*ones(N/b,1));
    rind = reshape(rind,[1, M*N/b^2]);
    rind = reshape(cell2mat(rind), [b^2, M*N/b^2]);
else 
    T = R;
    tind = 1:numel(T);
    rind = tind;
end

end

function Y = bin2(X)
% simple 2-fold binning
[m, n] = size(X);

Y = X( (2:2:m)-1,   (2:2:n)-1 ) +...
    X( (2:2:m)  ,   (2:2:n)-1 ) +...
    X( (2:2:m)-1,   (2:2:n)   ) +...
    X( (2:2:m)  ,   (2:2:n)   );

end

