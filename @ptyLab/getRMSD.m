function getRMSD(obj, positionIndex)
% getRMSD:      root mean square deviation between ptychogram and intensity estimate (I)

obj.params.currentDetectorError = abs(obj.params.Imeasured - obj.params.Iestimated);

if obj.params.saveMemory
    
    if obj.params.FourierMaskSwitch && ~obj.params.CPSCswitch
        obj.params.errorAtPos(positionIndex) = sum2(obj.params.currentDetectorError .* obj.params.W);
    elseif obj.params.FourierMaskSwitch && obj.params.CPSCswitch
         obj.params.errorAtPos(positionIndex) = ...
             sum( abs( sum( obj.params.Iestimated(         obj.params.upsampledIndex), 1 ) - ...
                           obj.params.ImeasuredDownsampled(obj.params.downsampledIndex) ) .* obj.params.W(:)' );
    else
        obj.params.errorAtPos(positionIndex) = sum2(obj.params.currentDetectorError);
    end

else
    
    obj.params.detectorError(:,:,positionIndex) = obj.params.currentDetectorError;

end