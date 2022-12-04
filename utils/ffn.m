function G = ffn(g)
% G = ffn(g)
% ffn performs n-dimensional Fourier transformation
% ffn is normalized (i.e. norm(g) = norm(G) )

G = fftshift( fftn( ifftshift( g ) ) ) / sqrt(numel( g ));

end