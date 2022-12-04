function [R, C] = GenerateArchimedeanSpiral(n, ds, varargin)
% [r, c] = GenerateArchimedeanSpiral(n, ds)
% optional: 
% [R, C] = GenerateArchimedeanSpiral(n, ds, varargin)
% ex.: [R, C] = GenerateArchimedeanSpiral(5, 4, 'a', 3, 'b', 1, 'c', 1);
% n  ~number of cycles of spiral
% ds ~number of pixels (converted to arc length) between each scan point
% R: row-position vector
% C: column-position vector
% a, b, c: parameter to control form of Archimedean spiral:
% r = b*theta.^(1/c)
% see also: https://en.wikipedia.org/wiki/Archimedean_spiral

% parse optional inputs
p = inputParser;
p.addParameter('b', 5)
p.addParameter('c', 1)
p.parse(varargin{:})

% generate angular coordinate
dtheta = ds / (2 * pi * n);
theta = 0 : dtheta : 2*n*pi;
% theta = sqrt(theta);

% for k = 1:n
%     dtheta = ds/(2*pi*k);
%     theta = [theta, (2*(k-1)*pi : dtheta : 2*k*pi)];
% end

% gemerate radial coordinate
r = p.Results.b * theta.^(1/p.Results.c);

% generate column (= x-direction) and row (= y-direction) position vectors
% C = round(r .* cos(theta));
% R = round(r .* sin(theta));
C = (r .* cos(theta));
R = (r .* sin(theta));

% for k = 2:length(C)
%     d = (C(k) - C(k-1))
%     if sqrt(  )
% end

R = [0,R,0];
C = [0,C,0];

end