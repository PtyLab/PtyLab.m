function g = ifft2c(G)
% G = ifft2c(g)
% ifft2c performs 2-dimensional inverse Fourier transformation
% ifft2c is normalized (i.e. norm(g) = norm(G) ), 
% i.e. it preserves the L2-norm
% if G is two-dimensional, ifft2c(G) yields the 2D iDFT of G
% if G is tree-dimensional, ifft2c(G) yields the 2D iDFT of G for each slice
% along the third dimension

% g = fftshift(ifft2(ifftshift(G))) * sqrt(numel(G));
g = fftshift(ifft2(ifftshift(G))) * sqrt( numel(G(:,:,1)) );
% note that numel is only executed along a single slice, so that it is
% independent of the number of slices! (important for correct normalization 
% under partially coherent conditions)
end
