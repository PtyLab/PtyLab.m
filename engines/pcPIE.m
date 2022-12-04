function obj = pcPIE(obj)
% mPIE: ePIE with accelerated gradient

obj.mPIEoperations('mode', 'init')

% mPIE parameter
t = 0;          % counter
beta = 0.3;     % feedback
eta = 0.7;      % friction (0.9)

% set object and probe buffers
obj.params.objectBuffer = obj.object;
obj.params.probeBuffer  = obj.probe;

% position correction parameters
daleth = 250;      % feedback (notice this is called alpha in manuscript)
beth = 0.9;        % friction
d = zeros(obj.numFrames, 2); % position search direction

% predefine shifts
rowShifts = [-1,-1,-1, 0, 0, 0, 1, 1, 1]';
colShifts = [-1, 0, 1,-1, 0, 1,-1, 0, 1]';

startAtIteration = 20;

% % weighting matrix
% W = ones(obj.Np, obj.Np, 'like', obj.probe);
% W(1,:) =  0;
% W(end,:)= 0;
% W(:,1) =  0;
% W(:,end)= 0;

for loop = 1:obj.params.numIterations
    
    % set position order
    obj.params.positionIndices = setPositionOrder(obj);
    
    % probe power correction
    if obj.params.probePowerCorrectionSwitch
        obj.probe = obj.probe / sqrt( sum( obj.probe(:) .* conj(obj.probe(:)) ) ) * obj.params.probePowerCorrection;
    end
    
    for positionLoop = 1:obj.numFrames
        
        % time increment
        t = t + 1;
        
        % get current position index
        positionIndex = obj.params.positionIndices(positionLoop);
        
        % get object patch
        row = obj.positions(positionIndex,1);
        col = obj.positions(positionIndex,2);
        objectPatch = obj.object(row:row+obj.Np-1, col:col+obj.Np-1,:);
        % note that object patch has size of probe array
        
        % form exit surface wave
        obj.params.esw = bsxfun(@times, objectPatch, obj.probe );
        
        % intensity constraint
        intensityProjection(obj, positionIndex);
        
        % difference term
        DELTA = obj.params.eswUpdate - obj.params.esw;
        
        % regularization
        aleph = 0.1;
         
        % object update
        absP2 = abs(obj.probe).^2;
        if strcmp(obj.operationMode, 'FPM')
            Pmax = max(max(sum(absP2, 3)));
             frac = abs(obj.probe) ./ Pmax .* conj(obj.probe) ./ ( aleph * max(max(sum(absP2,3))) + (1-aleph) * absP2 );
        else
            frac = conj(obj.probe) ./ ( aleph * max(max(sum(absP2,3))) + (1-aleph) * absP2 );
        end
        obj.object(row:row+obj.Np-1, col:col+obj.Np-1, :) = objectPatch + obj.params.betaObject * sum( bsxfun(@times, DELTA, frac), 3);
        
        % probe update
        absO2 = abs(objectPatch).^2;
        Omax = max(abs(obj.object(:)));
        frac = conj(objectPatch) ./ ( aleph * Omax + (1-aleph) * absO2 );
        obj.probe = obj.probe + obj.params.betaProbe * bsxfun(@times, DELTA, frac );
        
        % position correction
        if length(obj.params.error) > startAtIteration
            % position gradients
            % shift function (should/can this be preallocated?)
            shiftFcn = @(r,c) circshift(objectPatch, [r c]);
            shiftedImages = arrayfun(shiftFcn, rowShifts, colShifts, 'UniformOutput', false);
            shiftedImages = cat(3, shiftedImages{:});
            
            % truncated cross-correlation
            cc = squeeze( sum2( bsxfun(@times, conj(shiftedImages), ...
                obj.object(row:row+obj.Np-1, col:col+obj.Np-1)) ));
            cc = abs(cc);
            
            
            normFactor = real(sum2( conj(objectPatch) .* objectPatch ));
            grad_x = daleth * sum( (cc - mean(cc))/normFactor .* colShifts );
            grad_y = daleth * sum( (cc - mean(cc))/normFactor .* rowShifts );
            
            r = 3;
            if abs(grad_x) > r
                grad_x = r * grad_x/abs(grad_x);
            end
            if abs(grad_y) > r
                grad_y = r * grad_y/abs(grad_y);
            end
            
            d(positionIndex,:) = gather([grad_y, grad_x]) + beth * d(positionIndex,:);
        end
        
        % momentum update
        if rand > 0.95
            
            % object update
            objectGradient = obj.params.objectBuffer - obj.object;
            obj.params.objectMomentum =  objectGradient + eta * obj.params.objectMomentum;
            obj.object = obj.object - beta * obj.params.objectMomentum;
            
            % probe update
            probeGradient = obj.params.probeBuffer - obj.probe;
            obj.params.probeMomentum = probeGradient + eta * obj.params.probeMomentum;
            obj.probe = obj.probe - beta * obj.params.probeMomentum;
            
            % get object buffer
            obj.params.objectBuffer = obj.object;
            obj.params.probeBuffer = obj.probe;
            
        end
        
    end
    
    if length(obj.params.error) > startAtIteration
        % update positions
        obj.positions = obj.positions - round( d );
        
        
        % fix center of mass of positions
        obj.positions(:,1) = obj.positions(:,1) - round( mean(obj.positions(:,1)) - mean(obj.positions0(:,1)) );
        obj.positions(:,2) = obj.positions(:,2) - round( mean(obj.positions(:,2)) - mean(obj.positions0(:,2)) );
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        figure(102)
        scatter(obj.positions0(:,1), obj.positions0(:,2),'bo','filled') 
        axis square, hold on
        scatter(obj.positions(:,1), obj.positions(:,2),'yo','filled')
        hold off, legend('correct', 'estimated')
        set(gca, 'FontSize', 20)
        title('estimated positions')
        
        figure(103)
        scatter(d(:,1), d(:,2),'mo','filled'), 
        axis([-r r -r r]), axis square
        title('displacement')
        xlabel('[px]'), ylabel('[px]')
        set(gca, 'FontSize', 20)
    end
    
    % get error metrics
    obj.getErrorMetrics
    
    % modulus enforced probe
    if obj.params.modulusEnforcedProbeSwitch
        % propagate probe to detector
        obj.params.esw = obj.probe;
        obj.object2detector
        
        if obj.params.FourierMaskSwitch
            obj.params.ESW = ...
                bsxfun(@times, obj.params.ESW .* sqrt(obj.params.emptyBeam ./ (1e-10+sum( abs(obj.params.ESW).^2 , 3 ) )), obj.params.W) + ...
                bsxfun(@times, obj.params.ESW, 1-obj.params.W);
        else
            obj.params.ESW = obj.params.ESW .* sqrt(obj.params.emptyBeam ./ (1e-10+sum( abs(obj.params.ESW).^2 , 3 ) ) );
        end
        
        obj.detector2object
        
        obj.probe = obj.params.esw;
    end
    
    % orthogonalize modes
    if mod(loop, obj.params.orthogonalizationFrequency) == 0
        if obj.params.npsm > 1
            obj.mPIEoperations('mode', 'orthogonalize')
        end
        
        % center of mass stabilization
        if obj.params.comStabilizationSwitch
            centerOfMassStabilization(obj);
        end
    end
    
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
    
    if obj.params.absorbingProbeBoundary
        if strcmp(obj.operationMode, 'CPM')
            aleph = 5e-2;
        else
            aleph = 100e-2;
        end
        obj.probe = (1-aleph) * obj.probe + aleph * bsxfun(@times, obj.probe, obj.params.probeWindow);
    end
    
    if obj.params.absObjectSwitch
        obj.object = (1-obj.params.absObjectBeta) * obj.object + ...
            obj.params.absObjectBeta * abs(obj.object);
    end
    
    if obj.params.objectContrastSwitch
        % this is intended to slowly push non-measured object region to 
        % abs value lower than the max abs inside object roiallowing for good contrast when
        % monitoring object
        obj.object = 0.995 * obj.object + 0.005 * mean(mean(abs(obj.object(obj.params.objectROI(1):obj.params.objectROI(2), obj.params.objectROI(3):obj.params.objectROI(4),1 )))); 
    end
    
    if mod(loop, obj.monitor.figureUpdateFrequency) == 0
        % show reconstruction
        showReconstruction(obj)
    end
    
end

obj.mPIEoperations('mode', 'clearMemory')

end