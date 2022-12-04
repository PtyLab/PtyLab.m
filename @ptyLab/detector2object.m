function detector2object(obj)
% propagate from detector to object plane

switch obj.propagator.type
        
    case 'Fraunhofer'
        obj.ifft2s;
        
    case 'Fresnel'
        obj.ifft2s;
        obj.params.esw = bsxfun(@times, obj.params.esw, conj(obj.propagator.quadraticPhase));
        obj.params.eswUpdate = bsxfun(@times, obj.params.eswUpdate, conj(obj.propagator.quadraticPhase));
    
    case {'aspw','ASP'}
        obj.params.eswUpdate = ifft2c(bsxfun(@times, fft2c(obj.params.ESW), conj(obj.propagator.transferFunction)));
        
    case 'scaledASP'
        obj.params.eswUpdate = bsxfun(@times, ifft2c(  bsxfun(@times, fft2c( obj.params.ESW ), conj(obj.propagator.Q2) ) ), conj(obj.propagator.Q1));
    
    otherwise
        error('obj.propagator.type not properly specified')
        
end

return