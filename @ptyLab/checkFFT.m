function checkFFT(obj)

if obj.params.fftshiftSwitch
    %
    if obj.params.fftshiftFlag == 0
        disp('check fftshift...')
        disp('fftshift data for fast far-field update')
        % shift detector quantities
        obj.ptychogram = ifftshift( ifftshift( obj.ptychogram, 1 ), 2);
        if isfield(obj.params, 'ptychogramDownsampled')
            obj.params.ptychogramDownsampled = ifftshift( ifftshift( obj.params.ptychogramDownsampled, 1 ), 2);
        end
        if isfield(obj.params, 'W')
            obj.params.W = ifftshift( ifftshift( obj.params.W, 1 ), 2);
        end
        if isfield(obj.params, 'empyBeam')
            obj.params.empyBeam = ifftshift( ifftshift( obj.params.empyBeam, 1 ), 2);
        end
        if isfield(obj.params,'PSD')
            obj.params.PSD = ifftshift( ifftshift( obj.params.PSD, 1 ), 2);
        end
        % set glag to 1 since now the data is shifted
        obj.params.fftshiftFlag = 1;
    end
else
    if obj.params.fftshiftFlag == 1
        disp('check fftshift...')
        disp('ifftshift data')
        obj.ptychogram = fftshift( fftshift( obj.ptychogram, 1 ), 2);
        if isfield(obj.params, 'ptychogramDownsampled')
            obj.params.ptychogramDownsampled = fftshift( fftshift( obj.params.ptychogramDownsampled, 1 ), 2);
        end
        if isfield(obj.params,'W')
            obj.params.W = fftshift( fftshift( obj.params.W, 1 ), 2);
        end
        if isfield(obj.params, 'empyBeam')
            obj.params.empyBeam = fftshift( fftshift( obj.params.empyBeam, 1 ), 2);
        end
        % set glag to 0 since now the data is not shifted
        obj.params.fftshiftFlag = 0;
    end
end