function obj = pSD(obj)
% parallel steepest descent (pSD) algorithm

% define index function to access all object patches (submatrices) at once
indexFcn = @(r,c) obj.object(r:(r+obj.Np-1), c:(c+obj.Np-1));

for loop = 1:obj.params.numIterations
    
    % get object patches
    objectPatches = arrayfun(indexFcn, obj.positions(:,1), obj.positions(:,2), 'UniformOutput', false);
    
    % reshape submatrices along 4rd dimension and multiply with probe (faster way than this??)
    objectPatches = cat(4, objectPatches{:});
    %     objectPatches = reshape(objectPatches{:}, [obj.Np,obj.Np,1,obj.numFrames]);
    
    % multiply with probe
    obj.params.esw = bsxfun( @times, objectPatches, obj.probe );
    %     obj.params.esw = objectPatches .* obj.probe;
    
    % propagate to detector plane
    obj.object2detector;
    
    % save intensity at random position for inspection
    obj.params.currentPosition = randi(obj.numFrames);
    obj.params.Iestimated = abs(obj.params.ESW(:,:,obj.params.currentPosition)).^2;
    obj.params.Imeasured = obj.ptychogram(:,:,obj.params.currentPosition);
    
    % standard intensity update (todo: extend to multiple spatial modes)
    %     obj.params.ESW = obj.params.ESW .* sqrt( obj.ptychogram ./ ( abs(obj.params.ESW).^2 + 1e-10) );
    obj.params.ESW = bsxfun(@times, obj.params.ESW, reshape(sqrt( obj.ptychogram ./ ( squeeze(sum(abs(obj.params.ESW).^2, 3)) + 1e-10) ), [obj.Np,obj.Np,1,obj.numFrames]));
    
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
    
    for k = 1:obj.numFrames
        row = obj.positions(k,1);
        col = obj.positions(k,2);
        
        objectUpdateNumerator(row:row+obj.Np-1,col:col+obj.Np-1,:) = objectUpdateNumerator(row:row+obj.Np-1,col:col+obj.Np-1,:) + abs(obj.probe) / max(abs(obj.probe(:))) .* conj(obj.probe) .* DELTA(:,:,:,k);
        objectUpdateDenominator(row:row+obj.Np-1,col:col+obj.Np-1) = objectUpdateDenominator(row:row+obj.Np-1,col:col+obj.Np-1,:) + real(sum(abs(obj.probe).^2,3));
        
        %         objectUpdateNumerator(row:row+obj.Np-1,col:col+obj.Np-1,:) = objectUpdateNumerator(row:row+obj.Np-1,col:col+obj.Np-1,:) + sum(conj(obj.probe) .* DELTA(:,:,:,k),3);
%         objectUpdateNumerator(row:row+obj.Np-1,col:col+obj.Np-1,:) = objectUpdateNumerator(row:row+obj.Np-1,col:col+obj.Np-1,:) + conj(obj.probe) .* DELTA(:,:,:,k);
%         objectUpdateDenominator(row:row+obj.Np-1,col:col+obj.Np-1) = objectUpdateDenominator(row:row+obj.Np-1,col:col+obj.Np-1,:) + real(sum(abs(obj.probe).^2,3));
    end
    
    % regularization for SD updates
    gimmel = 0.5;
    
    % probe update
    %     obj.probe = obj.probe +  obj.params.betaProbe * probeUpdateNumerator ./ (gimmel * max(probeUpdateDenominator(:)) + (1-gimmel) * probeUpdateDenominator );
    obj.probe = obj.probe + obj.params.betaProbe * bsxfun(@times, probeUpdateNumerator, 1 ./ (gimmel * max(probeUpdateDenominator(:)) + (1-gimmel) * probeUpdateDenominator ));
    
    
    if ~isreal(objectUpdateDenominator)
        error('complex found'),
    end
    obj.object = obj.object + obj.params.betaObject * sum( bsxfun(@times, objectUpdateNumerator, 1./  (gimmel * max(objectUpdateDenominator(:)) + (1-gimmel)*objectUpdateDenominator)) , 3 );
    
    % get error metrics
    obj.params.errorAtPos = squeeze( sum(sum2( abs(DELTA).^2 ),3) );
    obj.params.errorAtPos = obj.params.errorAtPos ./ (obj.params.energyAtPos + 1);
    eAverage = sum( obj.params.errorAtPos );
    
    % append to error vector (for ploting error as function of iteration)
    obj.params.error = [obj.params.error eAverage];
    
    if obj.params.objectSmoothenessSwitch
        % object smootheness regularization (only appleid inside relevant object ROI)
        h1 = gaussian2D(2*obj.params.objectSmoothenessWidth + 1, obj.params.objectSmoothenessWidth);
        obj.object(obj.params.objectROI(1):obj.params.objectROI(2), obj.params.objectROI(3):obj.params.objectROI(4)) = ...
            (...
            (1-obj.params.objectSmoothnessAleph) * abs(obj.object(obj.params.objectROI(1):obj.params.objectROI(2), obj.params.objectROI(3):obj.params.objectROI(4))) + ...
            obj.params.objectSmoothnessAleph * normconv2( abs( obj.object(obj.params.objectROI(1):obj.params.objectROI(2), obj.params.objectROI(3):obj.params.objectROI(4))) , h1)) .* ...
            exp(1i * angle(obj.object(obj.params.objectROI(1):obj.params.objectROI(2), obj.params.objectROI(3):obj.params.objectROI(4))));
    end
    
    if obj.params.probeSmoothenessSwitch
        % probe regularization
        h2 = gaussian2D(2*obj.params.probeSmoothenessWidth + 1, obj.params.probeSmoothenessWidth);
        for k = 1:obj.params.npsm
            obj.probe(:,:,k) = (1-obj.params.probeSmoothnessAleph) * obj.probe(:,:,k) + ...
                obj.params.probeSmoothnessAleph * normconv2(obj.probe(:,:,k), h2);
        end
    end
    
    %     if obj.params.probeTVregSwitch
    %         epsilon = 1e-2;
    %         Gr = grad( obj.probe );
    %         d = sum3(Gr.^2,3);
    %         G = -div( Gr ./ repmat( sqrt( epsilon^2 + d ) , [1 1 2]) );
    %         obj.probe = obj.probe - 1e-2 * G;
    %     end
    
    if obj.params.absorbingProbeBoundary
        aleph = 100e-2;
        obj.probe = (1-aleph) * obj.probe + aleph * bsxfun(@times, obj.probe, obj.params.probeWindow);
    end
    
    % orthogonalize (spatial) modes
    if mod(loop, obj.params.orthogonalizationFrequency) == 0
        obj.orthogonalize
        
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