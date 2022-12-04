function [R, C] = GenerateConcentricGrid(Nr, s, rend)
% generate concentric circles
% Nr: number of circles (or shells)
% s: number of pixels between
% rend: end radius size (in pixel units)

dx = 1;     % max Resolution (Schritt von einem zum anderen Pixel)
rstart = dx;

r = linspace(rstart,rend,Nr);

% determine number of positions on k'th shell
nop = zeros(Nr,1);
for k = 1 : Nr
    nop(k) = floor( 2*pi*r(k) / s );
end        

positions = zeros(sum(nop), 2);
ind = 1;
for k = 1 : Nr
    
    dtheta = 2 * pi / nop(k);
    theta = (1:nop(k)) * dtheta + 2*pi/k;
    
    for l = 1:nop(k)
        positions(ind,:) = r(k) * [cos(theta(l)), sin(theta(l))];
        ind = ind + 1;
    end
end

positions = floor(positions / dx);% + 200 * ones(size(positions));

R = positions(:,1)';
C = positions(:,2)';
R = [0, R];
C = [0, C];

end