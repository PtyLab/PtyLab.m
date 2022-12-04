function getErrorMetrics(obj)
% get various error metrics

if ~obj.params.saveMemory
    
    % calculate mean error for all positions (make separate function for all of that)
    if obj.params.FourierMaskSwitch
        obj.params.errorAtPos = squeeze(sum(sum( ...
            bsxfun(@times, abs(obj.params.detectorError), obj.params.W), 1 ), 2));
    else
        obj.params.errorAtPos = squeeze(sum(sum( abs(obj.params.detectorError), 1 ), 2));
    end
    
end

obj.params.errorAtPos = obj.params.errorAtPos ./ (obj.params.energyAtPos + 1);
eAverage = sum( obj.params.errorAtPos );

% append to error vector (for ploting error as function of iteration)
obj.params.error = [obj.params.error eAverage];

end