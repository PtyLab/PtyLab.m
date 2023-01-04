function r = max_center( P )
% center quadratic probe function (or mode stack)

N = size(P,1);
x = linspace(-N/2,N/2-1,N);
% x = linspace(-N/2,N/2,N);
[X,Y] = meshgrid(x);

% form intensity
I = sum( abs(P).^2, 3 );
% normalize
I = I / sum2(I);

[row, col] = find( I == max(I(:)));
xs = X(row, col);
ys = Y(row, col);

Fx = X/N;
Fy = Y/N;

r = ifft2c(  bsxfun(@times, fft2c(P), exp(1i * 2*pi * (Fx * xs + Fy * ys) ) ) );







