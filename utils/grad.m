function del = grad( f )
% compute gradient of function using finite difference

delx = f( :, [2:end, end] ) - f;
dely = f( [2:end, end], : ) - f;
del = cat(3, delx, dely);

end
