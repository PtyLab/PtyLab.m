function obj = e3PIE(obj)
% multislice version of ePIE

if ~isfield(obj.params,'objectSlices')
    obj.params.objectSlices = bsxfun(@times, ones(obj.No, obj.No, obj.params.numSlices,'like',obj.object), nthroot(abs(obj.object), obj.params.numSlices).*exp(1i*angle(obj.object)/obj.params.numSlices) );
    
end

% initial temporal slices
probeTemp = bsxfun(@times, ones(obj.Np,obj.Np,obj.params.npsm, obj.params.numSlices,'like',obj.probe), obj.probe);
eswTemp = ones(obj.Np,obj.Np,obj.params.npsm, obj.params.numSlices,'like',obj.probe);

% preallocate transfer function
[~, H] = aspw( probeTemp(:,:,1,1), obj.params.dz, obj.lambda, obj.Lp);
% shift transfer function to avoid fftshifts for FFTs
H = ifftshift( ifftshift( H, 1 ), 2);


for loop = 1:obj.params.numIterations
    
    % set position order
    obj.params.positionIndices = setPositionOrder(obj);
    
    for positionLoop = 1:obj.numFrames
        
        positionIndex = obj.params.positionIndices(positionLoop);
        
        % get row and column indices of current position
        row = obj.positions(positionIndex,1);
        col = obj.positions(positionIndex,2);
        
        % form first exit surface wave
        probeTemp(:,:,:,1) = obj.probe;
        eswTemp(:,:,:,1) = bsxfun(@times, probeTemp(:,:,:,1), ...
            obj.params.objectSlices(row:row+obj.Np-1, col:col+obj.Np-1,1));
        
        % propagate through object
        for sliceLoop = 2:obj.params.numSlices
            probeTemp(:,:,:,sliceLoop) = ifft2( bsxfun(@times, fft2(eswTemp(:,:,:,sliceLoop-1)), H) );
            eswTemp(:,:,:,sliceLoop)  = bsxfun(@times, probeTemp(:,:,:,sliceLoop), ...
                obj.params.objectSlices(row:row+obj.Np-1, col:col+obj.Np-1, sliceLoop));
        end
        
%       % set esw behind last slice for FFT
        obj.params.esw = eswTemp(:,:,:,end);

        % intensity constraint (computes eswUpdate and other quantities)
        intensityProjection(obj, positionIndex);
        
        % calculate first difference term
        DELTA = obj.params.eswUpdate - obj.params.esw;
        
        % update object slices
        for loopTemp = 1:obj.params.numSlices-1
            sliceLoop = obj.params.numSlices - loopTemp + 1;
            % compute current object slice
            objectUpdate = objectPatchUpdate( obj,...
                obj.params.objectSlices(row:row+obj.Np-1, col:col+obj.Np-1, sliceLoop), DELTA, probeTemp(:,:,:,sliceLoop) );
            
            % eswTemp (here probe incident on last slice) update
            beth = 1;
            probeTemp(:,:,:,sliceLoop) = ...
                probeUpdate( obj, obj.params.objectSlices(row:row+obj.Np-1, col:col+obj.Np-1, sliceLoop), DELTA, probeTemp(:,:,:,sliceLoop), beth );
            
            % update object
            obj.params.objectSlices(row:row+obj.Np-1, col:col+obj.Np-1, sliceLoop) = objectUpdate;
            
            % back-propagate and calulate gradient term
            DELTA = ifft2( bsxfun(@times, fft2(probeTemp(:,:,:,sliceLoop)), conj(H) )) - eswTemp(:,:,:,sliceLoop-1);
        end
        
        % update last object slice
        objectUpdate = objectPatchUpdate( obj,...
                obj.params.objectSlices(row:row+obj.Np-1, col:col+obj.Np-1, 1), DELTA, probeTemp(:,:,:,1) );
            
        % update probe
        probeTemp(:,:,:,1) = probeUpdate( obj, obj.params.objectSlices(row:row+obj.Np-1, col:col+obj.Np-1, 1), DELTA, probeTemp(:,:,:,1), obj.params.betaProbe );
        
        % update last object
        obj.params.objectSlices(row:row+obj.Np-1, col:col+obj.Np-1, 1) = objectUpdate;
        
        % update probe
        obj.probe = probeTemp(:,:,:,1);
        
    end
    % set updated first object slice as product of all slices
    obj.object = prod(obj.params.objectSlices, 3);

    % get error metrics
    obj.getErrorMetrics
    
    % orthogonalize modes
%     if mod(loop, obj.params.orthogonalizationFrequency) == 0
%         error('mixmedode')
%         obj.orthogonalize
%     end
    
    % probe power correction
    if obj.params.probePowerCorrectionSwitch
        obj.probe = obj.probe / sqrt( sum( obj.probe(:) .* conj(obj.probe(:)) ) ) * obj.params.probePowerCorrection;
    end
    
    % object smootheness regularization
    if obj.params.objectSmoothenessSwitch
        % regularization
        h1 = fspecial('disk', obj.params.objectSmoothenessWidth);
%         gimmel = 5e-2;
        gimmel = obj.params.objectSmoothnessAleph;
        for k = 1:obj.params.numSlices
            obj.params.objectSlices(:,:,k) = (1-obj.params.objectSmoothnessAleph) * obj.params.objectSlices(:,:,k) + obj.params.objectSmoothnessAleph * convolve2(obj.params.objectSlices(:,:,k), h1, 'same');
            obj.params.objectSlices(:,:,k) = (1-gimmel) * obj.params.objectSlices(:,:,k) + gimmel * mean(mean(abs(obj.params.objectSlices(:,:,k)))) .* exp(1i * angle(obj.params.objectSlices(:,:,k)));
        end
    end
    
%     idx = find(abs(obj.params.objectSlices)>1);
%     obj.params.objectSlices(idx) = exp(1i * angle(obj.params.objectSlices(idx)));
    
    % center of mass stabilization
    if obj.params.comStabilizationSwitch
        centerOfMassStabilization(obj);
    end
    
    % show reconstruction
    if mod(loop, obj.monitor.figureUpdateFrequency) == 0 
        
        figure(89)
        switch obj.params.objectPlot
            case 'angle'
                
                imagesc(angle(obj.params.objectSlices(obj.params.objectROI(1):obj.params.objectROI(2), obj.params.objectROI(3):obj.params.objectROI(4),1))), colormap bone, axis image 
                subplot(1,2,2)
                imagesc(angle(obj.params.objectSlices(obj.params.objectROI(1):obj.params.objectROI(2), obj.params.objectROI(3):obj.params.objectROI(4),end))), colormap bone, axis image
            case 'abs'
                subplot(1,2,1)
                imagesc(abs(obj.params.objectSlices(obj.params.objectROI(1):obj.params.objectROI(2), obj.params.objectROI(3):obj.params.objectROI(4),1)),[0,1.2]), colormap bone, axis image off
                subplot(1,2,2)
                imagesc(abs(obj.params.objectSlices(obj.params.objectROI(1):obj.params.objectROI(2), obj.params.objectROI(3):obj.params.objectROI(4),end)),[0,1.2]), colormap bone, axis image off
                
            otherwise
                subplot(1,2,1)
                hsvxplot(obj.params.objectSlices(obj.params.objectROI(1):obj.params.objectROI(2), obj.params.objectROI(3):obj.params.objectROI(4),1), 'intensityScale', [0 1]) 
                subplot(1,2,2)
                hsvxplot(obj.params.objectSlices(obj.params.objectROI(1):obj.params.objectROI(2), obj.params.objectROI(3):obj.params.objectROI(4),end), 'intensityScale',[0 1])
        end
        subplot(1,2,1)
        title('slice 1'), set(gca, 'FontSize', 20), zoom(obj.params.objectZoom)
        subplot(1,2,2)
        title(['slice ',num2str(size(obj.params.objectSlices,3))]), set(gca, 'FontSize', 20), zoom(obj.params.objectZoom)
        showReconstruction(obj);
    end
    
end

end

%%% local functions %%%
function objectPatch = objectPatchUpdate( obj, objectPatch, DELTA, localProbe)


%     frac = conj(obj.probe) ./ max( max( sum( abs(obj.probe).^2, 3 ) ) );
frac = conj(localProbe) ./ max( max( sum( abs(localProbe).^2, 3 ) ) );
% see Thibault and Menzel 2013, nature, supplement for a derivation

n = obj.params.npsm;
%         n = max(1, obj.params.npsm-1);
objectPatch = objectPatch + obj.params.betaObject * sum( frac(:,:,1:n) .* DELTA(:,:,1:n), 3);
% CHANGE THIS OPTION OVER SWITCH
%         objectPatch = objectPatch + obj.params.betaObject * frac(:,:,1) .* DELTA(:,:,1);
%         objectPatch = objectPatch + betaObject * frac(:,:,1) .* DELTA(:,:,1);

end

function r = probeUpdate(obj, objectPatch, DELTA, localProbe, beth)

frac = conj(objectPatch) ./ max( max( sum( abs(objectPatch).^2, 3) ) );
if obj.params.npsm == 1
    % r = obj.probe + obj.params.betaProbe * sum( frac .* DELTA, 3 );
    % CHANGE THIS OPTION OVER SWITCH
    r = localProbe + beth * frac(:,:,1) .* DELTA(:,:,1);
else
    r = localProbe + beth * bsxfun(@times, DELTA, frac );
end

% 
% if obj.params.absorbingProbeBoundary
%     aleph = 1e-3;
%     r = (1-aleph) * r + aleph * bsxfun(@times, r, obj.params.probeWindow);
% end
end