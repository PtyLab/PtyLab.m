function obj = lsqPIE(obj)
% automatic step size ePIE
% as_ePIE
% see 
% Odstrčil, Michal, Andreas Menzel, and Manuel Guizar-Sicairos. 
% "Iterative least-squares solver for generalized maximum-likelihood ptychography." 
% Optics express 26.3 (2018): 3108-3123.
% note: the algorithm used here (lsqPIE) uses the adapted step size
% as reported in the above reference (for only a single scan position, not batches)
obj.params.stepProbeHistory = zeros(obj.params.npsm,obj.params.numIterations);
obj.params.stepObjectHistory = zeros(obj.params.numIterations,1);
for loop = 1:obj.params.numIterations
    
    % set position order
    obj.params.positionIndices = setPositionOrder(obj);
    
    stepObjectAverage = 0;
    stepProbeAverage = 0;
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
        chi = obj.params.eswUpdate - obj.params.esw; % chi
        
        % object update
        objectGrad = sum( conj(obj.probe) .* chi, 3 );
        
        % probe update
        probeGrad = conj(objectPatch) .* chi;
        
        % compute step size via full matrix equation
%         A = diag( [sum(sum2( abs(objectGrad .* obj.probe).^2 )); squeeze( sum2( abs(probeGrad .* objectPatch).^2 ))]  );
%         A(2:end,1) = squeeze(sum2(conj(probeGrad.*objectPatch).*obj.probe .* objectGrad));
%         A(1,2:end) = A(2:end,1)';
%         A = A+eye(size(A));
%         
%         b = [sum( sum2( real( conj(objectGrad .* obj.probe) .* chi  ) ) );...
%            squeeze(sum2( real( conj(probeGrad .* objectPatch) .* chi ) ))];
%         c = A\b;
%         stepSizeObject = c(1);
%         stepSizeProbe = reshape(c(2:end), [1,1,obj.params.npsm]);
        
        % diagonal matrix approximation 
        % compute step size probe
        gimmel = 0;
        stepSizeProbe = sum2( real( conj(probeGrad .* objectPatch) .* chi ) ) ./ ...
                       (sum2( abs(probeGrad .* objectPatch).^2 ) + gimmel);
        
        stepProbeAverage = stepProbeAverage + stepSizeProbe;
        
        % compute step size object
        daleth = 0;
        stepSizeObject = sum( sum2( real( conj(objectGrad .* obj.probe) .* chi  ) ) ) / ...
                        (sum(sum2( abs(objectGrad .* obj.probe).^2 )) + daleth);
        
        stepObjectAverage = stepObjectAverage + stepSizeObject;
        
        % probe update
        if strcmp(obj.operationMode, 'FPM')
            % FPM update
            obj.probe = obj.probe + stepSizeProbe .* abs(obj.probe) / max(abs(obj.probe(:))) .* probeGrad;
        else
            % CPM update
            obj.probe = obj.probe + stepSizeProbe .* probeGrad;
        end
        
        % set updated object patch
        obj.object(row:row+obj.Np-1, col:col+obj.Np-1,:) = ...
            objectPatch + stepSizeObject * objectGrad;
        
        
    end
    obj.params.stepObjectHistory(loop) = stepObjectAverage / obj.numFrames;
    obj.params.stepProbeHistory(:,loop) = squeeze(stepProbeAverage) / obj.numFrames;
    
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
            obj.params.objectSmoothnessAleph * normconv2( abs( obj.object(obj.params.objectROI(1):obj.params.objectROI(2), obj.params.objectROI(3):obj.params.objectROI(4))) , h1)) .* ...
            exp(1i * angle(obj.object(obj.params.objectROI(1):obj.params.objectROI(2), obj.params.objectROI(3):obj.params.objectROI(4))));
        
        % Good's roughness on probe
%         for k = 1:obj.params.npsm
%             
%             epsilon = 1e-2;
%             temp = obj.probe(:,:,k);
%             Gr = grad( temp );
%             G = -div( Gr ./ repmat( abs(temp) + epsilon , [1 1 2]) );
%             obj.probe(:,:,k) = obj.probe(:,:,k) - 20e-2 * G;
%           
%         end
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
    % r = obj.probe + obj.params.betaProbe * sum( frac .* DELTA, 3 );
    % CHANGE THIS OPTION OVER SWITCH
    r = obj.probe + obj.params.betaProbe * frac(:,:,1) .* DELTA(:,:,1);
else
    r = obj.probe + obj.params.betaProbe * bsxfun(@times, DELTA, frac );
end


if obj.params.absorbingProbeBoundary
    aleph = 1e-3;
    r = (1-aleph) * r + aleph * bsxfun(@times, r, obj.params.probeWindow);
end
end