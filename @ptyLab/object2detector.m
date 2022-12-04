function object2detector(obj)
% propagate from object to detector plane

switch obj.propagator.type
    
    case 'Fraunhofer'
        obj.fft2s;
        
    case 'Fresnel'
        obj.params.esw = bsxfun(@times, obj.params.esw, obj.propagator.quadraticPhase);
        obj.fft2s;
        
    case {'aspw','ASP'}
        obj.params.ESW = ifft2c(bsxfun(@times, fft2c(obj.params.esw), obj.propagator.transferFunction));
        
    case 'scaledASP'
        obj.params.ESW = ifft2c( bsxfun(@times, fft2c( bsxfun(@times, obj.params.esw, obj.propagator.Q1)  ), obj.propagator.Q2) );
        
    otherwise
        error('obj.propagator.type not properly specified')
        
end

return