function [Z, j] = zernike(r,theta,m,n)
% compute zernike polynomial with indices m, n
% m: azimuthal
% n: radial
% references:
% https://en.wikipedia.org/wiki/Zernike_polynomials
% difference in notation:
% notice the normalization by 1/sqrt(pi), which causes the set of all
% zernike polynomials as defined here to be orthonormal
if m == 0      % m zero
    Z = sqrt((n+1)/pi)*rfun(m,n,r);
else
    if m > 0 % m positive
        Z = sqrt(2*(n+1)/pi)*rfun(m,n,r) .* cos(m*theta);
    else     % m negative
        Z = sqrt(2*(n+1)/pi)*rfun(m,n,r) .* sin(m*theta);
    end
end
Z = Z .* (r<=1);
j = getOSAindex(m,n);
end

function R = rfun(m, n, r)
% radial factor of Zernike  polynomial
R = 0;
if mod(m-n,2) == 0
    for s = 0 : ((n-m)/2)
        R = R + r.^(n-2*s) * ...
            ( (-1)^s * gamma(n-s+1) ) / ...
            ( gamma(s+1) * gamma((n+m)/2-s+1) * gamma((n-m)/2-s+1) );
    end
else
    R = 0;
end

end

function j = getOSAindex(m,n)
% OSA index convention
% https://en.wikipedia.org/wiki/Zernike_polynomials
j = (n*(n+2) + m)/2;
end