function r = im2rgb( varargin )

switch nargin
    
    case 2
        im1 = varargin{1};
        im2 = varargin{2};
        
        % normalize
        im1 = round(im1/max(im1(:)));
        im2 = round(im2/max(im2(:)));
        
        r = uint16( cat(3, im1, im2, zeros(size(im1))) * 2^16 - 1);
end