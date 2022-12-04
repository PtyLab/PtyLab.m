function obj = PIE(obj)
% original PIE engine (with probe recovery)

for loop = 1:obj.params.numIterations
    
    % set position order
    obj.params.positionIndices = setPositionOrder(obj);
        
    for positionLoop = 1:obj.numFrames
        
        % get current position index
        positionIndex = obj.params.positionIndices(positionLoop);
        
        % get object patch
        row = obj.positions(positionIndex,1);
        col = obj.positions(positionIndex,2);
        objectPatch = obj.object(row:row+obj.Np-1, col:col+obj.Np-1,:);
        % note that object patch has size of probe array
        
        % form exit surface wave
        obj.params.esw = bsxfun(@times, objectPatch, obj.probe );
        
        % intensity constraint (computes eswUpdate and other quantities)
        intensityProjection(obj, positionIndex);
        
        % PIE
        % difference term
        DELTA = obj.params.eswUpdate - obj.params.esw;
        
        % object update
        objectPatch = objectPatchUpdate( obj, objectPatch, DELTA );
        
        % probe update
        obj.probe = probeUpdate( obj, objectPatch, DELTA );
        
        % set updated object patch
        obj.object(row:row+obj.Np-1, col:col+obj.Np-1,:) = objectPatch;
        
    end
    
    % get error metrics
    obj.getErrorMetrics
    
    % orthogonalize modes
%     obj.orthogonalization
    
    % probe power correction
    if obj.params.probePowerCorrectionSwitch
        obj.probe = obj.probe / sqrt( sum( obj.probe(:) .* conj(obj.probe(:)) ) ) * obj.params.probePowerCorrection;
    end
    
    % probe smootheness
    if obj.params.probeSmoothenessSwitch
        % probe regularization
        h2 = gaussian2D(2*obj.params.probeSmoothenessWidth + 1, obj.params.probeSmoothenessWidth);
        for k = 1:obj.params.npsm
            obj.probe(:,:,k) = (1-obj.params.probeSmoothnessAleph) * obj.probe(:,:,k) + ...
                obj.params.probeSmoothnessAleph * convolve2(obj.probe(:,:,k), h2, 'same');
        end
    end
    
    % TV
    if obj.params.probeTVregSwitch
        epsilon = 1e-2;
        for k = 1:obj.params.npsm;
            Gr = grad( obj.probe(:,:,k) );
            d = sum3(Gr.^2,3);
            G = -div( Gr ./ repmat( sqrt( epsilon^2 + d ) , [1 1 2]) );
            obj.probe(:,:,k) = obj.probe(:,:,k) - 1e-4 * G;
        end
    end
    
    if obj.params.comStabilizationSwitch
        % center of mass stabilization of probe
        centerOfMassStabilization(obj);
    end
    
    % show reconstruction
    if mod(loop, obj.monitor.figureUpdateFrequency) == 0, showReconstruction(obj); end
end

end


%%% local functions %%%

function objectPatch = objectPatchUpdate( obj, objectPatch, DELTA )

aleph = 1e-8;
absP = abs(obj.probe);
maxAbsP = max(max(sum(absP,3)));
frac = absP/maxAbsP .* conj(obj.probe) ./ (absP.^2 + aleph);
% see Thibault and Menzel 2013, nature, supplement for a derivation

objectPatch = objectPatch + obj.params.betaObject * sum(bsxfun(@times, DELTA, frac),3);

end


function r = probeUpdate(obj, objectPatch, DELTA)

aleph = 1e-8;
absO = abs(objectPatch);
maxAbsO = max(max(absO));
frac = absO/maxAbsO .* conj(objectPatch) ./ (absO.^2 + aleph);

r = obj.probe + obj.params.betaProbe * bsxfun(@times, DELTA, frac);

if obj.params.absorbingProbeBoundary
    r = bsxfun(@times, r, obj.params.probeWindow);
end
end