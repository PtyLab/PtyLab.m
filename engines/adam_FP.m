function obj = adam_FP(obj)
% momentum-accelerated gradient descent

% define index function to access all object patches (submatrices) at once
indexFcn = @(r,c) obj.object(r:(r+obj.Np-1), c:(c+obj.Np-1));

obj.mPIEoperations('mode', 'init')

t = 0;
daleth = 0.1;        % feedback
beth = 0.99;         % friction
s_probe = 1/2;             % adaptive step size
s_object = 1/2;
% set object and probe buffers
obj.params.objectBuffer = obj.object;
obj.params.probeBuffer  = obj.probe;

for loop = 1:obj.params.numIterations
    t = t+1;
    idx = randperm(obj.numFrames);
    idx = idx(1:obj.params.batchSize);
    
    % get object patches
    objectPatches = arrayfun(indexFcn, obj.positions(idx,1), obj.positions(idx,2), 'UniformOutput', false);
    
    % reshape submatrices along 4rd dimension and multiply with probe (faster way than this??)
    objectPatches = cat(4, objectPatches{:});
    %     objectPatches = reshape(objectPatches{:}, [obj.Np,obj.Np,1,obj.numFrames]);
    
    % multiply with probe
    obj.params.esw = bsxfun( @times, objectPatches, obj.probe );
    %     obj.params.esw = objectPatches .* obj.probe;
    
    % propagate to detector plane
    obj.object2detector;
    
    % save intensity at random position for inspection (here only a single, random position is required)
    obj.params.currentPosition = randi(obj.params.batchSize);
    obj.params.Iestimated = abs(obj.params.ESW(:,:,obj.params.currentPosition)).^2;
    obj.params.Imeasured = obj.ptychogram(:,:,obj.params.currentPosition);
    
    % standard intensity update (todo: extend to multiple spatial modes)
    obj.params.ESW = bsxfun(@times, obj.params.ESW, reshape(sqrt( obj.ptychogram(:,:,idx) ./ ( squeeze(sum(abs(obj.params.ESW).^2, 3)) + 1e-10) ), [obj.Np,obj.Np,1,obj.params.batchSize]));
    
    % back-propagate to object
    obj.detector2object;
    
    % formate gradient term (note that if memory issues arise, DELTA cis not needed; eswUpdate can be overwritten)
    DELTA = obj.params.eswUpdate - obj.params.esw;
    
    % probe update
    %     probeUpdateNumerator = sum( DELTA .* conj(objectPatches) , 4);
    probeUpdateNumerator = sum( bsxfun(@times, DELTA, conj(objectPatches)) , 4);
    probeUpdateDenominator = real(sum( objectPatches .* conj( objectPatches ), 4 ));
    
    % preallocate large object gradient
    objectUpdateNumerator = zeros([size(obj.object),obj.params.npsm], 'like', obj.probe);
    objectUpdateDenominator = real(zeros(size(obj.object), 'like', obj.probe));
    
    for k = 1:length(idx)
        row = obj.positions(idx(k),1);
        col = obj.positions(idx(k),2);
        
        objectUpdateNumerator(row:row+obj.Np-1,col:col+obj.Np-1,:) = objectUpdateNumerator(row:row+obj.Np-1,col:col+obj.Np-1,:) + abs(obj.probe) / max(abs(obj.probe(:))) .* conj(obj.probe) .* DELTA(:,:,:,k);
        objectUpdateDenominator(row:row+obj.Np-1,col:col+obj.Np-1) = objectUpdateDenominator(row:row+obj.Np-1,col:col+obj.Np-1,:) + real(sum(abs(obj.probe).^2,3));
        
    end
    
    % regularization for SD updates
    gimmel = 0.01;
    
    % probe update
    obj.probe = obj.probe + obj.params.betaProbe * bsxfun(@times, probeUpdateNumerator, 1 ./ (gimmel * max(probeUpdateDenominator(:)) + (1-gimmel) * probeUpdateDenominator ));
    
    % object update 
    if ~isreal(objectUpdateDenominator)
        error('complex found'),
    end
    obj.object = obj.object + obj.params.betaObject * sum( bsxfun(@times, objectUpdateNumerator, 1./  (gimmel * max(objectUpdateDenominator(:)) + (1-gimmel)*objectUpdateDenominator)) , 3 );
    
    % apply momentum
    % object update
    if mod(t, 5) == 0
        
        % object update
        objectGradient = obj.params.objectBuffer - obj.object;
        obj.params.objectMomentum =  daleth*objectGradient + (1-daleth) * obj.params.objectMomentum;
        s_object = sqrt(beth * s_object + (1-beth) * norm(objectGradient,2)^2);
%         obj.object = obj.object - s * obj.params.objectMomentum;
        obj.object = obj.object - obj.params.betaObject/(s_object + eps) * obj.params.objectMomentum;
        
        % probe update
        probeGradient = obj.params.probeBuffer - obj.probe;
        obj.params.probeMomentum =  daleth*probeGradient + (1-daleth) * obj.params.probeMomentum;
        s_probe = sqrt(beth * s_probe + (1-beth) * norm(probeGradient(:),2)^2);
%         obj.probe = obj.probe - s * obj.params.probeMomentum;
        obj.probe = obj.probe - obj.params.betaProbe/(s_probe + eps) * obj.params.probeMomentum;
        
        % get object buffer
        obj.params.objectBuffer = obj.object;
        obj.params.probeBuffer = obj.probe;
        
    end
    
    % get error metrics
    obj.params.errorAtPos = obj.params.errorAtPos;
    obj.params.errorAtPos(idx) = squeeze( sum(sum2( abs(DELTA).^2 ),3) );
    obj.params.errorAtPos = obj.params.errorAtPos ./ (obj.params.energyAtPos + 1);
    eAverage = sum( obj.params.errorAtPos );
    
    % append to error vector (for ploting error as function of iteration)
    obj.params.error = [obj.params.error eAverage];
    
    if obj.params.objectSmoothenessSwitch
        % object smootheness regularization (only appleid inside relevant object ROI)
        h1 = gaussian2D(2*obj.params.objectSmoothenessWidth + 1, obj.params.objectSmoothenessWidth);
        temp = ifft2c(obj.object);
        temp = (1-obj.params.objectSmoothnessAleph) * temp + ...
               obj.params.objectSmoothnessAleph * convolve2(temp, h1, 'same');
        obj.object = fft2c(temp);
    end
    
    if obj.params.probeSmoothenessSwitch
        % probe regularization
        h2 = gaussian2D(2*obj.params.probeSmoothenessWidth + 1, obj.params.probeSmoothenessWidth);
        for k = 1:obj.params.npsm
            obj.probe(:,:,k) = (1-obj.params.probeSmoothnessAleph) * obj.probe(:,:,k) + ...
                obj.params.probeSmoothnessAleph * convolve2(obj.probe(:,:,k), h2, 'same');
        end
    end
    
    if obj.params.probeTVregSwitch
        epsilon = 1e-2;
        for k = 1:obj.params.npsm
            Gr = grad( obj.probe(:,:,k) );
            d = sum3(Gr.^2,3);
            G = -div( Gr ./ repmat( sqrt( epsilon^2 + d ) , [1 1 2]) );
            obj.probe(:,:,k) = obj.probe(:,:,k) - 1e-2/obj.numFrames * obj.params.batchSize * G;
        end
    end
    
    if obj.params.objectTVregSwitch
        epsilon = 1e-2;
        temp = ifft2c(obj.object);
        Gr = grad( temp );
        d = sum3(Gr.^2,3);
        G = -div( Gr ./ repmat( sqrt( epsilon^2 + d ) , [1 1 2]) );
        obj.object = obj.object - 1e-3 * fft2c(G);
    end
    
    if obj.params.absorbingProbeBoundary
        aleph = 100e-2;
        obj.probe = (1-aleph) * obj.probe + aleph * bsxfun(@times, obj.probe, obj.params.probeWindow);
    end
    
    % orthogonalize modes and center pupil
    if mod(loop, obj.params.orthogonalizationFrequency) == 0
        if obj.params.npsm > 1
            obj.mPIEoperations('mode', 'orthogonalize')
        end
        
        % center of mass stabilization
        if obj.params.comStabilizationSwitch
            centerOfMassStabilization(obj);
        end
    end
    
    % show reconstruction
    if mod(loop, obj.monitor.figureUpdateFrequency) == 0
        showReconstruction(obj)
    end
    
    
end