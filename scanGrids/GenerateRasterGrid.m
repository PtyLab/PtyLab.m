function [R, C] = GenerateRasterGrid(n, ds, varargin)
% function to generate raster grid
% [R, C] = GenerateRasterGrid(n, ds, varargin)
% n: number of points per dimension

% parse optional inputs
p = inputParser;
p.addParameter('type', 'randomOffset')
p.addParameter('amplitude', 1)
p.parse(varargin{:})

[I, J] = ind2sub(n * [1, 1], 1:n^2); % index to 2D coordinate conversion
C = (I * ds) - ds; % positions in physical units
R = (J * ds) - ds; % positions in physical units

if mod(n,2)==0
    C = (C - n * ds/2);
    R = (R - n * ds/2);
else
    C = (C - (n-1) * ds/2);
    R = (R - (n-1) * ds/2);
end

switch p.Results.type
    case 'randomOffset'
        C = C + round( p.Results.amplitude * ( -1 + 2 * rand(size(C))));
        R = R + round( p.Results.amplitude * ( -1 + 2 * rand(size(R))));
    case 'raster'
        % do nothing
    otherwise
        error('grid type not properly specified')
end

distance = R.^2 + C.^2;
[~, index] = min(distance);

% R(index(1)) = [];
% C(index(1)) = [];
% R = [0, R];
% C = [0, C];
R = round(R - mean(R));
C = round(C - mean(C));

end


