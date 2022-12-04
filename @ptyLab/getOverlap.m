function getOverlap(obj, ind1, ind2)
% get area overlap 
% getOverlap(obj, m, n) calculates overlap between m-th and n-th position 
% note: the area overlap is only calculated for main mode of probe
% 
% example: getOverlap(obj, 1, 2) calculates overlap between first and second scan position 

% calculate shifts (positions is row-column order, shifts are xy order)
sy = abs( obj.positions(ind2,1) - obj.positions(ind1,1) ) * obj.dxp;
sx = abs( obj.positions(ind2,2) - obj.positions(ind1,2) ) * obj.dxp;

% task 1: get linear overlap
obj.getBeamWidth 
obj.params.linearOverlap = 1 - sqrt(single(sx^2 + sy^2)) / (min( obj.params.beamWidthX, obj.params.beamWidthY ));
obj.params.linearOverlap = max( obj.params.linearOverlap, 0 );
display(['estimated linear overlap: ',num2str(round(obj.params.linearOverlap * 1000)/10),' %'])

% task 2: get area overlap
% spatial frequency pixel size
df = 1/(obj.Np * obj.dxp);  
% spatial frequency meshgrid
[Fx, Fy] = meshgrid((-obj.Np/2:obj.Np/2-1)*df);
% absolute value of probe and 2D fft
p = abs(obj.probe(:,:,1));
PROBE = fft2c(p);
% calculate overlap between positions
obj.params.areaOverlap = abs( sum( sum( PROBE .* PROBE .* exp(-1i * 2*pi * ( Fx * sx + Fy * sy ) )  )) ) /...
                     sum( sum( abs(PROBE).^2 ));
display(['estimated area overlap: ',num2str(round(obj.params.areaOverlap * 1000)/10),' %'])
end