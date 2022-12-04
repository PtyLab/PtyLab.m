function r = normconv2(X, H)
% 2D convolution with energy conservation
% X: data to be 2D convolved
% H: convolution kernel

r = conv2(X, H, 'same');
r = r * norm(X,'fro') / norm(r, 'fro');
end


