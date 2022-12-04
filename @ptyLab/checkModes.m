function checkModes(obj)
% checkModes:       adjust probe and object array 
%dimensions if changed during reconstruction

% resize probe array
if size(obj.probe, 3) < obj.params.npsm
    % in this case, add probe modes
    
    addNum = obj.params.npsm - size(obj.probe, 3);
    
    disp('check modes...')
    disp('> probe mode structure changed')
    disp(['> add ',num2str(addNum),' probe modes'])
    disp(['> number of probe states: ',num2str(obj.params.npsm)])
    
    if obj.params.gpuFlag == 0
        obj.probe = cat(3, obj.probe, 1/10 * bsxfun(@times, rand(obj.Np, obj.Np, addNum, 'single'), obj.probe(:,:,1)));
    else
        obj.probe = cat(3, obj.probe, 1/10 * bsxfun(@times, gpuArray(rand(obj.Np, obj.Np, addNum, 'single')), obj.probe(:,:,1)));
    end
    probeModesChanged = true;
    
elseif size(obj.probe, 3) > obj.params.npsm
    % in this case, keep only the strongest probe modes
    
    discardNum = size(obj.probe, 3) - obj.params.npsm;
    
    disp('check modes...')
    disp('> probe mode structure changed')
    disp(['> discard ',num2str(discardNum),' probe modes'])
    disp(['> number of probe states: ',num2str(obj.params.npsm)])
    
    obj.probe = obj.probe(:,:,1:obj.params.npsm);
    
    probeModesChanged = true;
else
    % in this case, do not  change probe mode structure
    probeModesChanged = false;
end

if probeModesChanged
    % renormalize of probe energy
    obj.probe = obj.probe / sqrt( sum( obj.probe(:) .* conj(obj.probe(:)) ) ) * obj.params.probePowerCorrection;
end

% resize object array
if size(obj.object, 3) < obj.params.nosm
    
    addNum = obj.params.nosm - size(obj.object, 3);
    
    disp('check modes...')
    disp('> object mode structure changed')
    disp(['> add ',num2str(addNum),' object modes'])
    disp(['> number of object states: ',num2str(obj.params.nosm)])
    
    if obj.params.gpuFlag == 0
        obj.object = cat(3, obj.probe, 1/10*bsxfun(@times, rand(obj.No, obj.No, addNum, 'single'), obj.object));
    else
        obj.object = cat(3, obj.probe, 1/10*bsxfun(@times, gpuArray(rand(obj.No, obj.No, addNum, 'single')), obj.object));
    end
    
end

% resize 3D object array if required
if isfield(obj.params,'objectSlices')
    if strcmp(obj.params.engine,'e3PIE')
        if size(obj.params.objectSlices, 3) < obj.params.numSlices
            
            addNum = obj.params.numSlices - size(obj.params.objectSlices, 3);
            
            disp('check object slices...')
            disp('> object slice structure changed')
            disp(['> add ',num2str(addNum),' object slices'])
            disp(['> number of object states: ',num2str(obj.params.numSlices)])
            
            obj.params.objectSlices = cat(3, ones(obj.No, obj.No, addNum,'like',obj.object), obj.params.objectSlices);
            
        elseif size(obj.params.objectSlices, 3) > obj.params.numSlices
            
            discardNum = size(obj.params.objectSlices, 3) - obj.params.numSlices;
            
            disp('check object slices...')
            disp('> object slice structure changed')
            disp(['> discard ',num2str(discardNum),' object slices'])
            disp(['> number of object states: ',num2str(obj.params.numSlices)])
            
            obj.params.objectSlices = obj.params.objectSlices(:,:,1:obj.params.numSlices);
        end
        
    end
end

end