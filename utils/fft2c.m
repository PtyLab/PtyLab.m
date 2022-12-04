function G = fft2c(g)
% G = fft2c(g)
% fft2c performs 2-dimensional Fourier transformation
% fft2c is normalized (i.e. norm(g) = norm(G) ), 
% i.e. it preserves the L2-norm
% if g is two-dimensional, fft2c(g) yields the 2D iDFT of g
% if g is multi-dimensional, fft2c(g) yields the 2D iDFT of g for each slice
% along the third dimension


% G = fftshift(fft2(ifftshift(g))) / sqrt(numel(g));
G = fftshift(fft2(ifftshift(g))) / sqrt(numel( g(:,:,1)) );
% note that numel is only executed along a single slice, so that it is
% independent of the number of slices! (important for correct normalization 
% under partially coherent conditions)
end
