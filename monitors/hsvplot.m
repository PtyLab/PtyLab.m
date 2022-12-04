function varargout = hsvplot(varargin)
% hsvplot(u) generates hue-brightness plot of two dimensional input u 
% last change: 3rd March 2018

if nargin == 1
    u = varargin{1};
elseif nargin == 2
    x = varargin{1};
    u = varargin{2};
elseif nargin == 3
    x = varargin{1};
    y = varargin{2};
    u = varargin{3};
end

u = gather(u);

% normalize birghtness (value) to range [0,1]
r = abs(u);
r = r / ( max(r(:)) + eps );

% normalize angle 
phi = angle( u );
phi = ( phi + pi )/( 2 * pi );

% normalization of phase saturation
B = zeros(size(r, 1), size(r, 2), 3, 'like', r);         % Declare RGB array
B(:,:,1) = phi;
B(:,:,2) = 1;
B(:,:,3) = r;
A = hsv2rgb(B);

if nargin == 1
    imagesc(A); axis image 
elseif nargin == 2
    imagesc(x,x,A); axis image 
elseif nargin == 3
    imagesc(x,y,A); axis image 
end
set(gcf, 'Color', 'w');

switch nargout
    case 1
        varargout{1} = A;
    otherwise
end