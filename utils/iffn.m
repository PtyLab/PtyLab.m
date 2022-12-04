function g = iffn(G)
% G = ffn(g)
% iffn performs n-dimensional inverse Fourier transformation
% iffn is normalized (i.e. norm(g) = norm(G) )

g = fftshift( ifftn( ifftshift(G) ) ) * sqrt(numel(G));

end