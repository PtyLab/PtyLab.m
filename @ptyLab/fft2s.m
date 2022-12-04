function fft2s(obj)
% fft2c performs 2-dimensional Fourier transformation on 
% current exit wave
% fft2c is normalized (i.e. norm(g) = norm(G) ) 
% (for each slice if ESW is 3D stack) 
% i.e. it preserves the L2-norm (for each slice)
% if g is two-dimensional, fft2s(g) yields the 2D iDFT of g
% if g is multi-dimensional, fft2s(g) yields the 2D iDFT of 
% g for each slice along the third dimension

if obj.params.fftshiftSwitch
        obj.params.ESW = fft2(obj.params.esw) / obj.Np;
else 
    obj.params.ESW = fftshift(fft2(ifftshift(obj.params.esw))) / obj.Np;
end

% note that numel is only executed along a single slice, so that it is
% independent of the number of slices! (important for correct 
% normalization under partially coherent conditions)
end
