function [R, C] = GenerateNonUniformFermat(n, varargin)
% [R, C] = GenerateFermatSpiral(n, varargin)
% ex.: [R, C] = GenerateArchimedeanSpiral(n,'c', 1);
% R: row 
% C: column
% n: number of points generated
% c: optional argument that controls scaling of spiral
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