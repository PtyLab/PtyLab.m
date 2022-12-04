function checkGPU(obj, varargin)
% this function puts all relevant arrazs onto GPU

p = inputParser;
% transfer data only for specified position to GPU if required
p.addParameter('positionIndices', 1:obj.numFrames)
p.parse( varargin{:} );

if obj.params.gpuSwitch
    
    if obj.params.gpuFlag == 0
        disp('switch to gpu...')
        
        % clear gpu to prevent memory issues
        gpuDevice(1);
        
        % load properties to gpu
        obj.object = gpuArray(obj.object);
        obj.probe = gpuArray(obj.probe);
        if ~strcmp(obj.params.engine, 'kPIE')
            obj.ptychogram = gpuArray(obj.ptychogram);
        end
        obj.params.detectorError = gpuArray(obj.params.detectorError);
        obj.params.errorAtPos = gpuArray(obj.params.errorAtPos);
        obj.params.error = gpuArray(obj.params.error);
        
        % optional params
        if isfield(obj.params,'emptyBeam')
            obj.params.emptyBeam = gpuArray(obj.params.emptyBeam);
        end
        
        if isfield(obj.params,'background')
            obj.params.background = gpuArray(obj.params.background);
        end
        
        if isfield(obj.params,'probeStack')
            obj.params.probeStack = gpuArray( obj.params.probeStack );
        end
        
        if isfield(obj.params,'error')
            obj.params.error = gpuArray( obj.params.error );
        end
        
        % mPIE
        if isfield(obj.params,'objectMomentum')
            obj.params.objectMomentum = gpuArray( obj.params.objectMomentum );
        end
        
        if isfield(obj.params,'probeMomentum')
            obj.params.probeMomentum = gpuArray( obj.params.probeMomentum );
        end
        
        if isfield(obj.params,'objectBuffer')
            obj.params.objectBuffer = gpuArray( obj.params.objectBuffer );
        end
        
        if isfield(obj.params,'probeBuffer')
            obj.params.probeBuffer = gpuArray( obj.params.probeBuffer );
        end
        
        if isfield(obj.params,'objectSlices')
            obj.params.objectSlices = gpuArray( obj.params.objectSlices );
        end
        
        if isfield(obj.propagator,'quadraticPhase')
            obj.propagator.quadraticPhase = gpuArray( obj.propagator.quadraticPhase );
        end
        
        if isfield(obj.propagator,'transferFunction')
            obj.propagator.transferFunction = gpuArray( obj.propagator.transferFunction );
        end
        
        obj.params.gpuFlag = 1;
    end
    
else
    
    if obj.params.gpuFlag == 1
        disp('switch to cpu...')
        
        % load properties to cpu
        obj.object = gather(obj.object);
        obj.probe = gather(obj.probe);
        obj.ptychogram = gather(obj.ptychogram);
        obj.params.detectorError = gather(obj.params.detectorError);
        obj.params.errorAtPos = gather(obj.params.errorAtPos);
        obj.params.error = gather(obj.params.error);
        
        % optional params
        if isfield(obj.params,'emptyBeam')
            obj.params.emptyBeam = gather(obj.params.emptyBeam);
        end
        
        if isfield(obj.params,'background')
            obj.params.background = gather(obj.params.background);
        end
        
        if isfield(obj.params,'probeStack')
            obj.params.probeStack = gather( obj.params.probeStack );
        end
        
        if isfield(obj.params,'error')
            obj.params.error = gather( obj.params.error );
        end
        
        % mPIE
        if isfield(obj.params,'objectMomentum')
            obj.params.objectMomentum = gather( obj.params.objectMomentum );
        end
        
        if isfield(obj.params,'probeMomentum')
            obj.params.probeMomentum = gather( obj.params.probeMomentum );
        end
        
        if isfield(obj.params,'objectBuffer')
            obj.params.objectBuffer = gather( obj.params.objectBuffer );
        end
        
        if isfield(obj.params,'probeBuffer')
            obj.params.probeBuffer = gather( obj.params.probeBuffer );
        end
        
        if isfield(obj.params,'objectSlices')
            obj.params.objectSlices = gather( obj.params.objectSlices );
        end
        
        if isfield(obj.propagator,'quadraticPhase')
            obj.propagator.quadraticPhase = gather( obj.propagator.quadraticPhase );
        end
        
        if isfield(obj.propagator,'transferFunction')
            obj.propagator.transferFunction = gather( obj.propagator.transferFunction );
        end
        
        obj.params.gpuFlag = 0;
    end
    
end