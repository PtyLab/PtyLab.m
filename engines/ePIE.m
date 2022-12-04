function obj = ePIE(obj)
% extended ptychographic iterative engine

for loop = 1:obj.params.numIterations
    
    % set position order
    obj.params.positionIndices = setPositionOrder(obj);
    
    for positionLoop = 1:obj.numFrames
        
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
        
        % difference term
        DELTA = obj.params.eswUpdate - obj.params.esw;
        
        % object update
%         objectPatch = objectPatchUpdate( obj, objectPatch, DELTA );
        obj.object(row:row+obj.Np-1, col:col+obj.Np-1,:) = objectPatchUpdate( obj, objectPatch, DELTA );
        % probe update
        obj.probe = probeUpdate( obj, objectPatch, DELTA );
        
        % set updated object patch
%         obj.object(row:row+obj.Np-1, col:col+obj.Np-1,:) = objectPatch;
        
    end
    
    % get error metrics
    obj.getErrorMetrics
    
    % moduslus enforced probe
    if obj.params.modulusEnforcedProbeSwitch
        % propagate probe to detector
        obj.params.esw = obj.probe;
        obj.object2detector
        
        obj.params.ESW = ...
            bsxfun(@times, sqrt(obj.params.emptyBeam ./ (1e-10+sum( abs(obj.params.ESW).^2 , 3 ) )), obj.params.W) + ...
            bsxfun(@times, obj.params.ESW, 1-obj.params.W);
        
        obj.detector2object
        
        obj.probe = obj.params.esw;
    end
    
    % orthogonalize modes
    if mod(loop, obj.params.orthogonalizationFrequency) == 0
        obj.orthogonalize
    end
    
    % probe power correction
    if obj.params.probePowerCorrectionSwitch
        obj.probe = obj.probe / sqrt( sum( obj.probe(:) .* conj(obj.probe(:)) ) ) * obj.params.probePowerCorrection;
    end
    
    if obj.params.objectSmoothenessSwitch
        % object smootheness regularization (only appleid inside relevant object ROI)
        h1 = gaussian2D(2*obj.params.objectSmoothenessWidth + 1, obj.params.objectSmoothenessWidth);
        obj.object(obj.params.objectROI(1):obj.params.objectROI(2), obj.params.objectROI(3):obj.params.objectROI(4)) = ...
            (...
            (1-obj.params.objectSmoothnessAleph) * abs(obj.object(obj.params.objectROI(1):obj.params.objectROI(2), obj.params.objectROI(3):obj.params.objectROI(4))) + ...
            obj.params.objectSmoothnessAleph * convolve2( abs( obj.object(obj.params.objectROI(1):obj.params.objectROI(2), obj.params.objectROI(3):obj.params.objectROI(4))) , h1,'wrap')) .* ...
            exp(1i * angle(obj.object(obj.params.objectROI(1):obj.params.objectROI(2), obj.params.objectROI(3):obj.params.objectROI(4))));
        
    end
    
    if obj.params.probeSmoothenessSwitch
        % probe regularization
        h2 = gaussian2D(2*obj.params.probeSmoothenessWidth + 1, obj.params.probeSmoothenessWidth);
        for k = 1:obj.params.npsm
            obj.probe(:,:,k) = (1-obj.params.probeSmoothnessAleph) * obj.probe(:,:,k) + ...
                obj.params.probeSmoothnessAleph * convolve2(obj.probe(:,:,k), h2, 'same');
        end
    end
    
    if obj.params.absObjectSwitch
        obj.object = (1-obj.params.absObjectBeta) * obj.object + ...
                        obj.params.absObjectBeta * abs(obj.object);
    end
    
    if obj.params.comStabilizationSwitch
        % center of mass stabilization of probe
        centerOfMassStabilization(obj);
    end
    
    if obj.params.objectContrastSwitch
        % this is intended to slowly push non-measured object region to 
        % abs value lower than the max abs inside object roiallowing for good contrast when
        % monitoring object
        obj.object = 0.995 * obj.object + 0.005 * mean(mean(abs(obj.object(obj.params.objectROI(1):obj.params.objectROI(2), obj.params.objectROI(3):obj.params.objectROI(4),1 )))); 
    end
    
    % show reconstruction
    if mod(loop, obj.monitor.figureUpdateFrequency) == 0, showReconstruction(obj); end
    
end

end

%%% local functions %%%
function objectPatch = objectPatchUpdate( obj, objectPatch, DELTA)

frac = conj(obj.probe) ./ max( max( sum( abs(obj.probe).^2, 3 ) ) );
% see Thibault and Menzel 2013, nature, supplement for a derivation

if obj.params.nosm == 1
%     n = obj.params.npsm;
    %         n = max(1, obj.params.npsm-1);
    %         objectPatch = objectPatch + obj.params.betaObject * sum( frac(:,:,1:n) .* DELTA(:,:,1:n), 3);
    % CHANGE THIS OPTION OVER SWITCH
    objectPatch = objectPatch + obj.params.betaObject * frac(:,:,1) .* DELTA(:,:,1);
else
    objectPatch = objectPatch + obj.params.betaObject * bsxfun(@times, DELTA, frac);
end


end

function r = probeUpdate(obj, objectPatch, DELTA)


frac = conj(objectPatch) ./ max( max( sum( abs(objectPatch).^2, 3) ) );
if obj.params.npsm == 1
    r = obj.probe + obj.params.betaProbe * sum( frac .* DELTA, 3 );
    % CHANGE THIS OPTION OVER SWITCH
%     r = obj.probe + obj.params.betaProbe * frac(:,:,1) .* DELTA(:,:,1);
else
    r = obj.probe + obj.params.betaProbe * bsxfun(@times, DELTA, frac );
end


if obj.params.absorbingProbeBoundary
    aleph = 1e-3;
    r = (1-aleph) * r + aleph * bsxfun(@times, r, obj.params.probeWindow);
end
end