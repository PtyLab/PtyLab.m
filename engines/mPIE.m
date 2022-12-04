function obj = mPIE(obj)
% mPIE: ePIE with accelerated gradient
% see:  https://www.osapublishing.org/optica/abstract.cfm?uri=optica-4-7-736

obj.mPIEoperations('mode', 'init')

t = 0;          % counter
beta = 0.3;     % feedback
eta = 0.7;      % friction (0.9)

% set object and probe buffers
obj.params.objectBuffer = obj.object;
obj.params.probeBuffer  = obj.probe;

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
        
        if obj.params.objectTVregSwitch
            
            % TV regularization
            epsilon = 1e-2;
            temp = obj.object(row:row+obj.Np-1, col:col+obj.Np-1, :);
            Gr = grad( temp );
            d = sum(abs(Gr).^2,3);
            G = -div( Gr ./ repmat( sqrt( epsilon^2 + d ) , [1 1 2] ) );
            switch obj.operationMode
                case 'FPM'
                    error('todo')
                case 'CPM'
                    obj.object(row:row+obj.Np-1, col:col+obj.Np-1, :) = ...
                        obj.object(row:row+obj.Np-1, col:col+obj.Np-1, :) - obj.params.objectTVregStepSize * G;
            end
        end
        
        % probe update
        absO2 = abs(objectPatch).^2;
        Omax = max(abs(obj.object(:)));
        frac = conj(objectPatch) ./ ( aleph * Omax^2 + (1-aleph) * absO2 );
        obj.probe = obj.probe + obj.params.betaProbe * bsxfun(@times, DELTA, frac );
        
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
    
    % center of mass stabilization
    if obj.params.comStabilizationSwitch
        centerOfMassStabilization(obj);
    end
    
    % orthogonalize modes
    if mod(loop, obj.params.orthogonalizationFrequency) == 0
        if obj.params.npsm > 1
            obj.mPIEoperations('mode', 'orthogonalize')
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
        obj.object = (1-obj.params.absObjectStepSize) * obj.object + ...
            obj.params.absObjectStepSize * abs(obj.object);
    end
    
    if mod(loop, obj.monitor.figureUpdateFrequency) == 0
        % show reconstruction
        showReconstruction(obj)
    end
    
end

obj.mPIEoperations('mode', 'clearMemory')

end

