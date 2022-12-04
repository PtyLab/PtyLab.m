function Uout = two_step_prop(Uin, plane)
% two step propagator
% allows variable sampling between source and detector plane
% 
% (todo: check sampling conditions)

global lambda
global dxs zs dxp dxo zo dxd Xp Yp Xd Yd

switch plane
        case 'sourcePinhole'
        d1 = dxs;
        d2 = dxp;
        Dz = zs;
    case 'objectDetector_freeSpace'
        d1 = dxo;
        d2 = dxd;
        Dz = zo;
    otherwise
        error('propagation plane not specified')
end

N = size(Uin, 1);

% number of grid points
k = 2*pi/lambda; % wavenumber

% % source-plane coordinates
% [x1, y1] = meshgrid((-N/2 : 1 : N/2 - 1) * d1);

% magnification
m = d2/d1;

% intermediate plane
Dz1 = Dz / (1 - m); % propagation distance
d1a = lambda * abs(Dz1) / (N * d1);

% coordinates
[x1a, y1a] = meshgrid((-N/2 : N/2-1) * d1a);

% evaluate the Fresnel-Kirchhoff integral
QitmOut = 1 / (1i*lambda*Dz1) .* exp(1i*k/(2*Dz1) * (x1a.^2+y1a.^2));
QitmIn = exp(1i * k/(2*Dz1) * (Xp.^2 + Yp.^2));
Uitm =  bsxfun(@times, ft2(bsxfun(@times, Uin, QitmIn), d1), QitmOut);

% observation plane
Dz2 = Dz - Dz1; % propagation distance

% coordinates
% [x2, y2] = meshgrid((-N/2 : N/2-1) * d2);

% evaluate the Fresnel diffraction integral
Qout = 1 / (1i*lambda*Dz2) .* exp(1i*k/(2*Dz2) * (Xd.^2+Yd.^2));
Qin = exp(1i * k/(2*Dz2) * (x1a.^2 + y1a.^2));
% Uout = Qout .* ft2(Uitm .* Qin, d1a);
Uout =  bsxfun(@times, ft2(bsxfun(@times, Uitm, Qin), d1), Qout);


