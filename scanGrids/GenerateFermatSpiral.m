function [R, C] = GenerateFermatSpiral(n, varargin)
% [R, C] = GenerateFermatSpiral(n, varargin)
% ex.: [R, C] = GenerateArchimedeanSpiral(n,'c', 1);
% R: row 
% C: column
% n: number of points generated
% c: optional argument that controls scaling of spiral
% see https://en.wikipedia.org/wiki/Fermat%27s_spiral

% parse optional inputs
p = inputParser;
p.addParameter('c', 1) % constant multiplier
p.parse(varargin{:})

% golden ratio
r = p.Results.c * sqrt(1:n);
theta0 = 137.508/180*pi;
theta = (1:n) * theta0;

% generate angular coordinate
% theta = [];
% for k = 1:n
%     dtheta = ds/(2*pi*k);
%     theta = [theta, (2*(k-1)*pi : dtheta : 2*k*pi)];
% end

% gemerate radial coordinate
% r = sqrt( p.Results.c * theta );


% generate column (= x-direction) and row (= y-direction) position vectors
% C = round(r .* cos(theta));
% R = round(r .* sin(theta));
C = (r .* cos(theta));
R = (r .* sin(theta));

% R(end) = [];
% C(end) = [];
R = [0, R];
C = [0, C];

end