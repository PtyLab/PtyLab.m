function [R, C] = generate_non_uniform_fermat_spiral(n, varargin)
% [R, C] = generate_non_uniform_fermat_spiral(n, varargin)
% R: row 
% C: column
% n: number of points generated
% radius: optional argument that controls radius of spiral
% power: optional argument that controls non-uniform scaling of spiral
% see https://en.wikipedia.org/wiki/Fermat%27s_spiral

% parse optional inputs
p = inputParser;
p.addParameter('radius', 1000) % radius in micrometer
p.addParameter('power', 1)     % determines uniformness*
p.parse(varargin{:})

% * note:
% * p = 1 is standard Fermat 
% * p > 1 yields more points towards the center of grid  

% golden ratio
r = sqrt(1:n) / sqrt(n);
theta0 = 137.508 / 180 * pi;
theta = (1:n) * theta0;

C = p.Results.radius * r.^p.Results.power .* cos(theta);
R = p.Results.radius * r.^p.Results.power .* sin(theta);

R = [0, R];
C = [0, C];

end