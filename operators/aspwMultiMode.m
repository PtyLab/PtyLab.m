function u = aspwMultiMode(u, z, lambda, L)
%   ASPW wave propagation
%   u = aspw(u,z,lambda,L)
%   u: field distribution at z = 0 (u is assumed to be square, i.e. N x N)
%   z: propagation distance
%   lambda: wavelength
%   L: total size [m] of 
%   following: Matsushima et al., "Band-Limited Angular Spectrum Method for
%   Numerical Simulation of Free-Space
%   Propagation in Far and Near Fields", Optics Express, 2009

k = 2*pi/lambda;
N = size(u,1);
[Fx, Fy] = meshgrid((-N/2:N/2-1)/L);
f_max = L/(lambda*sqrt(L^2+4*z^2));    
W = rect(Fx/(2*f_max)).*rect(Fy/(2*f_max));  
H = exp( 1i * k * z * sqrt(1-(Fx*lambda).^2-(Fy*lambda).^2) );
U = fouriert2(u);
u = ifouriert2( bsxfun(@times, U, H .* W) );

end

% auxiliary functions

function G = fouriert2(g)
% 
% G=ft2(g,delta_x,delta_y)
% ft2 performs two-dimensional Fourier transformation

G = fftshift( fftn( ifftshift(g) ) ) / sqrt(numel(g));

end

function g = ifouriert2(G)
% 
% g=ift2(G,delta_fx,delta_fy)
% ift2 performs two-dimensional inverse Fourier transformation

g = fftshift( ifftn( ifftshift(G) ) ) * sqrt(numel(G));

end
