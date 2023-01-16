function r = normconv2(X, H)
% 2D convolution with energy conservation
% X: data to be 2D convolved
% H: convolution kernel

% r = conv2(X, H, 'same');
% r = r * norm(X,'fro') / norm(r, 'fro');

energyBefore = norm(X(:), 2);
r = X;
for k = 1:size(X, 3)
    r(:,:,k) = conv2(X(:,:,k), H, 'same');
end
energyAfter = norm(X(:), 2);
r = r / energyAfter * energyBefore;
end


