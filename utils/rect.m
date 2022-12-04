function y = rect( x )
% function y = rect( x )
% note that discontinuities have "sharp edges"
y = double( abs( x ) < 1/2 );
end