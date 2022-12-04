function varargout = hsvplot(u)
% hsvplot(u) generates hue-brightness plot of two dimensional input u 
% last change: 3rd March 2018

% normalize birghtness (value) to range [0,1]
r = abs(u);
r = r / ( max(r(:)) + eps );

% normalize angle 
phi = angle( u );
phi = ( phi + pi )/( 2 * pi );

% normalization of phase saturation

B = zeros(size(u,1),size(u,2),3);         %Declare RGB array
B(:,:,1) = phi;
B(:,:,2) = 1;
B(:,:,3) = r;
A = hsv2rgb(B);

imagesc(A); axis image 
set(gcf, 'Color', 'w');

switch nargout
    case 1
        varargout{1}=uint8(A*255);
    otherwise
end