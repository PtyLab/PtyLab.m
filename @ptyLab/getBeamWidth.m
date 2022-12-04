function getBeamWidth(obj)
% calculates beam width defined as FWHM

P2 = sum(abs(obj.probe).^2, 3);
P2 = P2 / sum(sum( P2 ));
xMean = sum(sum( obj.Xp .* P2 ));
yMean = sum(sum( obj.Yp .* P2 ));

xVariance = sum(sum( (obj.Xp - xMean).^2 .* P2 ));
yVariance = sum(sum( (obj.Yp - yMean).^2 .* P2 ));

% the constant c converts to FWHM 
% see e.g. https://en.wikipedia.org/wiki/Full_width_at_half_maximum
c = 2 * sqrt(2 * log(2));
obj.params.beamWidthX = gather( c * sqrt(xVariance) );
obj.params.beamWidthY = gather( c * sqrt(yVariance) );

end