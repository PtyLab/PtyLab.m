function obj = mqNewton(obj)
% quasi-Newton 2nd order

aperture = circ(obj.Xp, obj.Yp, obj.entrancePupilDiameter);
obj.mPIEoperations('mode', 'init')

t = 0;
daleth = 0.9;        % feedback
beth = 0.99;         % friction
s_probe = 1/2;             % adaptive step size
s_object = 1/2;
% set object and probe buffers
obj.params.objectBuffer = obj.object;
obj.params.probeBuffer  = obj.probe;

for loop = 1:obj.params.numIterations
    t = t +1 ;
    % set position order
    obj.params.positionIndices = setPositionOrder(obj);
    
    for positionLoop = 1:obj.numFrames
        
        positionIndex = obj.params.positionIndices(positionLoop);
        
        % get object patch
        row = obj.positions(positionIndex,1);
        col = obj.positions(positionIndex,2);
        objectPatch = obj.object(row:(row+obj.Np-1), col:(col+obj.Np-1),:);

        % form exit surface wave
        obj.params.esw = bsxfun(@times, objectPatch, obj.probe );
        
        % intensity constraint (computes eswUpdate and other quantities)
        intensityProjection(obj, positionIndex);
        
%         % propagate to detector plane
%         obj.object2detector;
%         
%         % save intensity at random position for inspection (here only a single, random position is required)
%         obj.params.currentPosition = positionIndex;
%         obj.params.Iestimated = abs(obj.params.ESW).^2;
%         obj.params.Imeasured = obj.ptychogram(:,:,positionIndex);
% 
%         % standard intensity update (todo: extend to multiple spatial modes)
%         obj.params.ESW = bsxfun(@times, obj.params.ESW, sqrt( obj.ptychogram(:,:,positionIndex) ./ (abs(obj.params.ESW).^2 + 1e-10) ));
% 
%         % back-propagate to object
%         obj.detector2object;
        
        % difference term
        DELTA = obj.params.eswUpdate - obj.params.esw;
        Omax = max(abs(obj.object(:)));
        
        obj.object(row:(row+obj.Np-1), col:(col+obj.Np-1),:) = objectPatch + obj.params.betaProbe * abs(obj.probe) .* conj(obj.probe) .* DELTA ./ ...
                       (max(abs(obj.probe(:))) .* (abs(obj.probe).^2 + 1));

        obj.probe = obj.probe +  obj.params.betaObject * abs(objectPatch) .* conj(objectPatch) .* DELTA ./ ...
                    (Omax .* (abs(objectPatch).^2 + 1 )) .* aperture;


    end
    
    % apply momentum
    % object update
%     if mod(t, 5) == 0
        
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
        
%     end
    
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