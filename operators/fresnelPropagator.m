function [r, dq, q, Qx, Qy] = fresnelPropagator( u, z, lambda, L)

k = 2*pi/lambda;
%% source coordinates
N = size(u,1);
dx = L / N;
x = (-N/2:N/2-1) * dx;
[X, Y] = meshgrid(x);

%% target coordinates
dq = lambda * z / L;
q = (-N/2:N/2-1) * dq;
[Qx, Qy] = meshgrid(q);

%% phase terms
Qin = exp(1i * k/(2*z) * (X.^2 + Y.^2));
% Qout = exp(1i * k/(2*z) * (Qx.^2 + Qy.^2));

% r = Qout .* fft2c(Qin .* u);
r =fft2c(Qin .* u);

end