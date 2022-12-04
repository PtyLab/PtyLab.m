function ifft2s(obj)
% ifft2c performs inverse 2-dimensional Fourier transformation on current
% detector wave
% ifft2c is normalized (i.e. norm(g) = norm(G) ) (for each slice if ESW is 3D stack) 
% i.e. it preserves the L2-norm (for each slice)
% if g is two-dimensional, ifft2s yields the 2D iDFT
% if g is multi-dimensional, fft2s yields the 2D iDFT for each slice of the input (along third dimension)

if obj.params.fftshiftSwitch
    obj.params.eswUpdate = ifft2(obj.params.ESW) * obj.Np;
else 
    obj.params.eswUpdate = fftshift(ifft2(ifftshift(obj.params.ESW))) * obj.Np;
end



% note that numel is only executed along a single slice, so that it is
% independent of the number of slices! (important for correct normalization 
% under partially coherent conditions)
end
