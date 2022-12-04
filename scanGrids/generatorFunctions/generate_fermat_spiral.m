function [R, C] = generate_fermat_spiral(n, varargin)
% [R, C] = generate_fermat_spiral(n, 'scaling', 1);
% genrates Fermat scan grid with n points
% R: row 
% C: column
% n: number of points generated
% scaling: optional argument that controls scaling of spiral
% see https://en.wikipedia.org/wiki/Fermat%27s_spiral

% parse optional inputs
p = inputParser;
p.addParameter('scaling', 1) % constant multiplier
p.parse(varargin{:})

% golden ratio
r = p.Results.scaling * sqrt(1:n);
theta0 = 137.508/180*pi;
theta = (1:n) * theta0;

% generate column (= x-direction) and row (= y-direction) position vectors
C = (r .* cos(theta));
R = (r .* sin(theta));

R = [0, R]';
C = [0, C]';

end