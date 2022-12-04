function z = circ(x, y, d)
% function z = circ(x, y, d)
% x: meshgrid (increasing along columns)
% y: meshgrid (increasing along rows)
% d: diameter

r = sqrt(x.^2 + y.^2);
z = single(r<d/2);
z(r==d/2) = 0;

end