function orthogonalize(obj)
% orthogonalize:    Orthogonalize beam and/or object transmission function
% example: obj.orthogonalize 


% SVD-based orthogonalization
% {

% temp = obj.probe;
if obj.params.npsm > 1
    % SVD
    [obj.probe, obj.params.normalizedEigenvaluesProbe, obj.params.MSPVprobe] = orthogonalizeModes(obj.probe, obj.params.npsm);

    % get probe purity (= spatial coherence measure)
    obj.params.purity = sqrt(sum( obj.params.normalizedEigenvaluesProbe.^2 ));
    
    % phase Synchronization 
%     for k = 1:obj.params.npsm
%         [obj.probe(:,:,k), c] = phaseSynchronization(temp(:,:,k), obj.probe(:,:,k));
%         obj.params.MSPVprobe(k,:) = obj.params.MSPVprobe(k,:) * c;
%     end
end



% object orthogonalization
% if obj.params.nosm > 1
%     [obj.object, obj.params.normalizedEigenvaluesObject, obj.params.MSPVobject] = orthogonalizeModes(obj.object, obj.params.nosm);
% end

end