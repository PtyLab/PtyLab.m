function r = gaussian2D( window_size, standard_deviation)
% creates 2D gaussian weights inside window 
% window_size: side length of window in units of pixels
% standard_deviation: standard deviation in units of pixels

if mod(window_size, 2) == 0
    n = window_size;
else
    n = (window_size - 1)/2;
end

x = linspace(-n,n,window_size);
y = x';

% use broadcasting (x:rows, y:cols)
r = exp( -( x.^2 + y.^2 )/(2*standard_deviation^2) );
r = r/sum(r(:));

end


