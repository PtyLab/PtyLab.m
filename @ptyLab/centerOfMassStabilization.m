function centerOfMassStabilization(obj)
% centers the probe

% calculate center of mass
P2 = sum(abs(obj.probe).^2, 3);
denom = sum(sum( P2 )) * obj.dxp;
xc = round(sum2( obj.Xp .* P2 ) / denom);
yc = round(sum2( obj.Yp .* P2 ) / denom);

% shift only if necessary
if xc^2 + yc^2 > 1
    % shift probe
    for k = 1:obj.params.npsm
        
        obj.probe(:,:,k) = circshift(obj.probe(:,:,k), -[yc, xc]);
        
        if strcmp(obj.params.engine, 'mPIE')
            obj.params.probeMomentum(:,:,k) = circshift(obj.params.probeMomentum(:,:,k), -[yc, xc]);
            obj.params.probeBuffer(:,:,k) = circshift(obj.params.probeBuffer(:,:,k), -[yc, xc]);
        end
    end
    
    % shift object
    for k = 1:obj.params.nosm
        
        obj.object(:,:,k) = circshift(obj.object(:,:,k), -[yc, xc]);
        
        if strcmp(obj.params.engine, 'mPIE')
            obj.params.objectMomentum(:,:,k) = circshift(obj.params.objectMomentum(:,:,k), -[yc, xc]);
            obj.params.objectBuffer(:,:,k) = circshift(obj.params.objectBuffer(:,:,k), -[yc, xc]);
        end
        
    end
end