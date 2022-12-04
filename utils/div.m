function r = div( del )

r = del(:,:,1)-del(:,[1 1:end-1],1) + ...
    del(:,:,2)-del([1 1:end-1],:,2);

end