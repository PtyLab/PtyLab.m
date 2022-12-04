function [Z,m,n] = zernikeNoll(r,theta,j)
% compute zernike polynomial with indices m, n
% m: azimuthal
% n: radial
% references:
% https://en.wikipedia.org/wiki/Zernike_polynomials
% difference in notation:
% notice the normalization by 1/sqrt(pi), which causes the set of all
% zernike polynomials as defined here to be orthonormal

[m,n] = convertj(j);

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

function [m,n] = convertj(k)
% OSA index convention
% https://en.wikipedia.org/wiki/Zernike_polynomials
% m: radial index
% n: azimuthal index

switch k
    case 1
        n = 0;
        m = 0;
    case 2
        n = 1;
        m = 1;
    case 3
        n = 1;
        m = -1;
    case 4
        n = 2;
        m = 0;
    case 5
        n = 2;
        m = -2;
    case 6
        n = 2;
        m = 2;
    case 7
        n = 3;
        m = -1;
    case 8
        n = 3;
        m = 1;
    case 9
        n = 3;
        m = -3;
    case 10
        n = 3;
        m = 3;
    case 11
        n = 4;
        m = 0;
    case 12
        n = 4;
        m = 2;
    case 13
        n = 4;
        m = -2;
    case 14
        n = 4;
        m = 4;
    case 15
        n = 4;
        m = -4;
    case 16
        n = 5;
        m = 1;
    case 17
        n = 5;
        m = -1;
    case 18
        n = 5;
        m = 3;
    case 19
        n = 5;
        m = -1;
    case 20
        n = 5;
        m = 5;
    case 21
        n = 5;
        m = -5;
    case 22
        n = 6;
        m = 0;
    case 23
        n = 6;
        m = -2;
    case 24
        n = 6;
        m = 2;
    case 25
        n = 6;
        m = -4;
end
end