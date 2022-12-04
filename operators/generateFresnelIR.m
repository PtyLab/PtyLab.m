function [hIRin, Qout] = generateFresnelIR(plane)
% 

global lambda
global Ns Np No Nl Nd
global dxs xs Xs Ys zs
global dxp xp Xp Yp zp
global dxo xo Xo Yo zo 
global dxl xl Xl Yl zl 
global dxd xd Xd Yd

k = 2*pi / lambda;

switch plane
    case 'sourcePinhole'
        hIRin = exp( 1i * k/(2*zs) * (Xs.^2 + Ys.^2) );
        Qout = exp( 1i * k/(2*zs) * (Xp.^2 + Yp.^2) );
    case 'objectDetector_freeSpace'
        hIRin = exp( 1i * k/(2*zo) * (Xo.^2 + Yo.^2) );
        Qout = exp( 1i * k/(2*zo) * (Xd.^2 + Yd.^2) );
    otherwise
        error('plane for Fresnel propagator not specified')
end



