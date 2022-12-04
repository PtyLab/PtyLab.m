function obj = zPIE(obj)
% zPIE

if ~isfield( obj.params, 'zHistory' )
    obj.params.zHistory = [];
end

zMomentum = 0;

% preallocate grids
if strcmp( obj.propagator.type, 'aspw' )
    n = obj.Np;
else
    n = 2*obj.Np;
end

[X, Y] = meshgrid((-n/2:n/2-1));w = exp( - (sqrt(X.^2 + Y.^2) / (obj.Np) ).^4 );

for loop = 1:obj.params.numIterations
    
    % set position order
    obj.params.positionIndices = setPositionOrder(obj);
    
    % get positions
    
    if loop == 1
        zNew = obj.zo;
    else
        d = 10;
        dz = linspace(-d * obj.params.DoF, d * obj.params.DoF, 11) / 10;
        merit = [ ];

        for k = 1:length(dz)
            imProp = aspw(w .* obj.object(end/2-n/2+1:end/2+n/2, end/2-n/2+1:end/2+n/2), dz(k), obj.wavelength, n*obj.dxo);
            
            % TV approach
            aleph = 1e-2;
            gradx = imProp - circshift(imProp, [0 1]);
            grady = imProp - circshift(imProp, [1 0]);
            merit = [merit, sum2( sqrt( abs(gradx).^2 + abs(grady).^2 + aleph ) )];    % c = 1 in paper
            
        end
        
        dz = sum(dz .* merit) / sum(merit);
        eta = 0.7;
        zMomentum = eta * zMomentum + obj.params.zPIEgradientStepSize * dz;
        zNew = obj.zo + zMomentum;
    end
    
    obj.params.zHistory = [obj.params.zHistory, gather(obj.zo)];
    
    disp(['position updated: z = ', num2str(obj.zo * 1e3), ' mm (dz = ', num2str(round(zMomentum * 1e7)/10), ' um)' ])
    
    % reset coordinates
    obj.zo = zNew;
    
    % resample
    if ~strcmp(obj.propagator.type, 'aspw')
        obj.dxo = obj.wavelength * obj.zo / obj.Ld;
        obj.positions = round( obj.params.encoder / obj.dxo );
        obj.positions = obj.positions + obj.No/2 - round(obj.Np/2);
        
        % object coordinates
%         obj.No = 2^12 + 2^11;
        obj.Lo = obj.No * obj.dxo;
        obj.xo = (-obj.No/2:obj.No/2-1) * obj.dxo;
        [obj.Xo, obj.Yo] = meshgrid(obj.xo);
        
        % probe coordinates
        obj.dxp = obj.dxo;
        obj.Np = obj.Nd;
        obj.Lp = obj.Np * obj.dxp;
        obj.xp = (-obj.Np/2:obj.Np/2-1)*obj.dxp;
        [obj.Xp, obj.Yp] = meshgrid(obj.xp);
        
        % reset propagator
        obj.propagator.quadraticPhase = exp(1i * pi/(obj.wavelength*obj.zo) * (obj.Xp.^2 + obj.Yp.^2));
    end
    
    for positionLoop = 1:obj.numFrames
        
        positionIndex = obj.params.positionIndices(positionLoop);
        
        %%% patch 1 %%%
        % get object patch1
        row1 = obj.positions(positionIndex,1);
        col1 = obj.positions(positionIndex,2);
        objectPatch1 = obj.object(row1:row1+obj.Np-1, col1:col1+obj.Np-1,:);
        
        % form exit surface wave1
        obj.params.esw = bsxfun(@times, objectPatch1, obj.probe );
        
        % intensity constraint (computes eswUpdate and other quantities)
        intensityProjection(obj, positionIndex);
        
        % difference term1
        DELTA = obj.params.eswUpdate - obj.params.esw;
        
        % object update
        objectPatch = objectPatchUpdate( obj, objectPatch1, DELTA );
        
        % set updated object patch
        obj.object(row1:row1+obj.Np-1, col1:col1+obj.Np-1,:) = objectPatch;
        
        % probe update
        obj.probe = probeUpdate( obj, objectPatch1, DELTA );
        
    end
    
    % get error metrics
    obj.getErrorMetrics
    
    % orthogonalize modes
    if mod(loop, obj.params.orthogonalizationFrequency) == 0
        obj.orthogonalize
    end
    
    % probe power correction
    if obj.params.probePowerCorrectionSwitch
        obj.probe = obj.probe / sqrt( sum( obj.probe(:) .* conj(obj.probe(:)) ) ) * obj.params.probePowerCorrection;
    end
    
    % object smootheness regularization
    if obj.params.objectSmoothenessSwitch
        % regularization
        h1 = fspecial('disk', obj.params.objectSmoothenessWidth);
        obj.object(obj.params.objectROI(1):obj.params.objectROI(2), obj.params.objectROI(3):obj.params.objectROI(4)) = ...
            (...
            (1-obj.params.objectSmoothnessAleph) * abs(obj.object(obj.params.objectROI(1):obj.params.objectROI(2), obj.params.objectROI(3):obj.params.objectROI(4))) + ...
            obj.params.objectSmoothnessAleph * convolve2( abs( obj.object(obj.params.objectROI(1):obj.params.objectROI(2), obj.params.objectROI(3):obj.params.objectROI(4))) , h1,'wrap')) .* ...
            exp(1i * angle(obj.object(obj.params.objectROI(1):obj.params.objectROI(2), obj.params.objectROI(3):obj.params.objectROI(4))));
        
        obj.object = 0.99 * obj.object + 0.01 * (mean(abs(obj.object(:))));
    end
    
    if obj.params.probeSmoothenessSwitch
        h2 = fspecial('disk', obj.params.probeSmoothenessWidth);
        for k = 1:obj.params.npsm
            obj.probe(:,:,k) = (1-obj.params.probeSmoothnessAleph) * obj.probe(:,:,k) + ...
                obj.params.probeSmoothnessAleph * convolve2(obj.probe(:,:,k), h2, 'same');
        end
    end
    
    if obj.params.absObjectSwitch
        obj.object = (1-obj.params.absObjectBeta) * obj.object + ...
            obj.params.absObjectBeta * abs(obj.object);
    end
    
    % center of mass stabilization
    if obj.params.comStabilizationSwitch
        centerOfMassStabilization(obj);
    end
    
    if obj.params.objectContrastSwitch
        % this is intended to slowly push non-measured object region to 
        % abs value lower than the max abs inside object roiallowing for good contrast when
        % monitoring object
        obj.object = 0.995 * obj.object + 0.005 * mean(mean(abs(obj.object(obj.params.objectROI(1):obj.params.objectROI(2), obj.params.objectROI(3):obj.params.objectROI(4),1 )))); 
    end
    
    % show reconstruction
    if mod(loop, obj.monitor.figureUpdateFrequency) == 0
        
        idx = linspace(0, log10(length(obj.params.zHistory)), min(length(obj.params.zHistory), 100));
        idx = round( 10.^idx);
        
        showReconstruction(obj);
        figure(666)
        set(gcf,'color','w')
        semilogx(idx, obj.params.zHistory(idx) * 1e3, 'o-')
        xlabel('iteration')
        ylabel('estimated z [mm]')
        set(gca,'Fontsize',20)
        axis square
        grid on
        drawnow
    end
    
end

end

%%% local functions %%%
function objectPatch = objectPatchUpdate( obj, objectPatch, DELTA)


frac = conj(obj.probe) ./ max( max( sum( abs(obj.probe).^2, 3 ) ) );
% see Thibault and Menzel 2013, nature, supplement for a derivation

if obj.params.nosm == 1
    %     n = obj.params.npsm;
            n = max(1, obj.params.npsm-1);
            objectPatch = objectPatch + obj.params.betaObject * sum( frac(:,:,1:n) .* DELTA(:,:,1:n), 3);
    % CHANGE THIS OPTION OVER SWITCH
%     objectPatch = objectPatch + obj.params.betaObject * frac(:,:,1) .* DELTA(:,:,1);
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
    aleph = 5e-2;
    r = (1-aleph) * r + aleph * bsxfun(@times, r, obj.params.probeWindow);
end
end