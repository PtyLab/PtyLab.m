function esw = inverseTwoStepPropagator(obj, ESW)

if size(ESW, 3) == 1
    % single mode implementation
    
    % invert phase term in detector plane
    ESW = ESW .* conj(obj.ODpropagator.QOut);
    
    % detector to intermediate plane 
    esw = ifft2c( ESW  ) .* conj( obj.ODpropagator.QIn );
    
    % invert phase term in intermediate plane
    esw = esw .* conj(obj.ODpropagator.QitmOut);
    
    % intermediate plane - object
    esw = ifft2c( esw ) .* conj( obj.ODpropagator.QitmIn );

%     esw = ifft2c(ESW);

elseif size(ESW, 3) > 1
    
    % multi mode implementation
    
    % invert phase term in detector plane
    ESW = bsxfun( @times, ESW, conj(obj.ODpropagator.QOut) );
    
    % detector to intermediate plane 
    esw = bsxfun(@times, ifft2c( ESW ), conj(obj.ODpropagator.QIn) );
    
    % invert phase term in intermediate plane
    esw = bsxfun( @times, esw, conj(obj.ODpropagator.QitmOut) );
    
    % intermediate plane - object
    esw = bsxfun(@times, ifft2c( esw ), conj(obj.ODpropagator.QitmIn) );

%     esw = ifft2c(ESW);
    
else
    error('this propagation case is not specified')
end


end