function intensityProjection(obj, positionIndex)
% intensity projection

% zero division mitigator
gimmel = 1e-10;

% propagate to detector
obj.object2detector;

% background mode
if and( obj.params.backgroundModeSwitch, not( isfield(obj.params, 'background') ) )
    obj.params.background = ones(obj.Nd, obj.Nd, 'like', obj.probe) * 1e-1;
end

% get estimated intensity
if obj.params.backgroundModeSwitch
    obj.params.Iestimated = real(sum( abs(obj.params.ESW).^2 , 3 ) + obj.params.background);
else
    obj.params.Iestimated = real(sum( abs(obj.params.ESW).^2 , 3 ));
end

% get measured intensity
if obj.params.CPSCswitch
    obj.decompressionProjection(positionIndex)
elseif and( strcmp(obj.params.engine, 'kPIE'), obj.params.gpuSwitch )
    obj.params.Imeasured = obj.params.diffractionData(:,:,positionIndex);
else
    obj.params.Imeasured = obj.ptychogram(:,:,positionIndex);
end

obj.params.currentPosition = positionIndex;

% todo: take out all PSDestimation!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% if obj.params.PSDestimationSwitch
%     % estimate PSD
%     obj.params.PSDestimate = obj.params.PSDestimate + obj.params.Iestimated;
% end

% calculate error
obj.getRMSD(positionIndex);
gimmel = 1;
lambda = 1;
% intensity projection
switch obj.params.intensityConstraint
    
    case 'Anscombe'
        frac = sqrt( (obj.params.Imeasured + gimmel) ./ (obj.params.Iestimated + gimmel) );
        
    case 'ProxAnscombe'
        frac = (sqrt(obj.params.Imeasured) + lambda * sqrt(obj.params.Iestimated) +  gimmel) ./ ((1+lambda) * sqrt(obj.params.Iestimated) + gimmel);
        
    case 'Poisson'
        frac = (obj.params.Imeasured + gimmel) ./ (obj.params.Iestimated + gimmel);
        
    case 'ProxPoisson'
        lambda = 1;
        frac = (obj.params.Imeasured + lambda * obj.params.Iestimated +  gimmel) ./ ((1+lambda) * obj.params.Iestimated + gimmel);
        
    case 'Gaussian'
        gimmel = 1;
        frac = 2 - (obj.params.Iestimated + gimmel) ./ (obj.params.Imeasured + gimmel);
    
    case 'fluctuation'
        
        % scaling
%         if obj.params.FourierMaskSwitch
%             aleph = sum(sum( obj.params.Imeasured .* obj.params.Iestimated .* obj.params.W )) / ...
%                     sum(sum( obj.params.Imeasured .* obj.params.Imeasured .* obj.params.W ));
%         else
%             aleph = sum(sum( obj.params.Imeasured .* obj.params.Iestimated )) / ...
%                     sum(sum( obj.params.Imeasured .* obj.params.Imeasured ));
%         end

        aleph = sum2(obj.params.Imeasured) / sum2(obj.params.Iestimated);
        
        obj.params.intensityScaling(positionIndex) = aleph;
        % scaled projection
%         frac = (1 + aleph)/2  * obj.params.Imeasured ./ (obj.params.Iestimated + gimmel);
%         frac = (1 + sqrt(aleph))/2  * obj.params.Imeasured ./ (obj.params.Iestimated + gimmel);
%         frac = sqrt((1 + aleph)/2 * obj.params.Imeasured ./ (obj.params.Iestimated + gimmel));
        daleth = 0.75;
        frac = sqrt( obj.params.Imeasured ./ ( ( 1-daleth + daleth * aleph) * obj.params.Iestimated + gimmel) );
        
    case 'exponential'
        
        x = obj.params.currentDetectorError ./ ...
            (obj.params.Iestimated + gimmel);
        % if the 1 is changed into larger values, does that mitigate background error?
        W = exp(-0.1 * x);
%         frac = sqrt( obj.params.Imeasured ./ (obj.params.Iestimated + gimmel) );
        frac = obj.params.Imeasured ./ (obj.params.Iestimated + gimmel);
        frac = W .* frac + (1-W);
        
    case 'MAP'
        frac = (obj.params.Imeasured+obj.params.Iestimated) ./ (2*obj.params.Iestimated + gimmel);
        
    otherwise
        
        % standard projection
%         frac = arrayfun(@sqrtFrac, obj.params.Imeasured, obj.params.Iestimated + gimmel);
        frac = sqrt(obj.params.Imeasured ./ (obj.params.Iestimated + gimmel));
        
        if ~strcmp(obj.params.intensityConstraint,'standard')
            warning('intensity constraint not properly specified!')
        end
end

% prevent numerical instability
frac(isnan(frac)) = 1;

% apply mask
if obj.params.FourierMaskSwitch && ~obj.params.CPSCswitch && length(obj.params.error) > 5
        frac = obj.params.W .* frac + (1-obj.params.W);
        

%     frac( ~obj.params.W ) = 1;
end

% update ESW
% if ~isreal(objectUpdateDenominator)
%     error('complex found'),
% end

if obj.params.gpuSwitch == 2
    obj.params.ESW = obj.params.ESW .* frac; % implicit expansion, compatible with CudaMat
else
    obj.params.ESW = bsxfun(@times, obj.params.ESW, frac ); % bsxfun: like implicit expansion but backwards compatible (< MATLAB2016b)
end

% update background
if obj.params.backgroundModeSwitch
    % see PhD thesis by Peng Li 
    if obj.params.FourierMaskSwitch
        obj.params.background = obj.params.background .* (1 + 1/obj.numFrames * (sqrt(frac) - 1) ).^2 .* obj.params.W;
    else
        obj.params.background = obj.params.background .* (1 + 1/obj.numFrames * (sqrt(frac) - 1) ).^2;
    end
end

% if obj.params.ESWSmoothenessSwitch
%     h1 = gaussian2D(2*obj.params.ESWSmoothenessWidth + 1, obj.params.ESWSmoothenessWidth);
%     for k = 1:obj.params.npsm
%         obj.params.ESW(:,:,k) = (1-obj.params.ESWSmoothnessAleph) * obj.params.ESW(:,:,k) + ...
%                             obj.params.ESWSmoothnessAleph * convolve2( obj.params.ESW(:,:,k), h1, 'wrap');
%     end
% end
% back propagate to object plane
obj.detector2object;

end