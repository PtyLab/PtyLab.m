function Uout = scaledASPinv(u, z, lambda, dx, dq)
% function [Uout] = ang_spec_prop(u, z, lambda, dx, dq)
% dx: grid spacing in original plane (u)
% dq: grid spacing in destination plane (Uout)
% 

% assume square grid
N = size(u,1);

% if (dq > 2*lambda*z/(N*dx))
%    error(' sampling grid in observation plane too coarse') 
% end

% optical wavenumber
k = 2*pi/lambda;

% source-plane coordinates
[x1, y1] = meshgrid((-N/2 : N/2 - 1) * dx);
r1sq = x1.^2 + y1.^2;
% spatial frequencies (of source plane)
df1 = 1 / (N*dx);
[fX, fY] = meshgrid((-N/2 : 1 : N/2 - 1) * df1);
fsq = fX.^2 + fY.^2;
% scaling parameter
m = dq/dx;

% quadratic phase factors
Q1 = exp(1i*k/2*(1-m)/z*r1sq);
Q2 = exp(-1i*pi^2*2*z/m/k*fsq);

% compute the propagated field
Uout = conj(Q1) .* ifft2c( conj(Q2) .* fft2c(  u ) );
