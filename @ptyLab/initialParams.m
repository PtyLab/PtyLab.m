function obj  = initialParams( obj, varargin )

% parse optional inputs
p = inputParser;
p.addParameter('absorbingProbeBoundary', false) % set probe to zero outside 2 x exit pupil diameter
p.addParameter('absObjectSwitch', false) % impose object being real and positiv
p.addParameter('absObjectStepSize', 0.05) % step size for absObjectSwitch
p.addParameter('batchSize', 10)         % fft-block size for parallel engines
p.addParameter('betaObject',0.25)       % mPIE feedback parameter (object)
p.addParameter('betaProbe',0.25)        % mPIE feedback parameter (probe)
p.addParameter('comStabilizationSwitch', true)% if true, the probe is centered
p.addParameter('deleteFrames',[])       % 1D vector with numbers indicating which frames to exclude
p.addParameter('engine','ePIE')         % plot every ... iterations
p.addParameter('fftshiftSwitch', true)  % switch to prevent excessive fftshifts
p.addParameter('figureUpdateFrequency',1) % plot every ... iterations
p.addParameter('FourierMaskSwitch',false) % switch to apply Fourier mask
p.addParameter('fontSize', 17)          % fontSize in plots
p.addParameter('initialObject','ones')  % choose initial object ('ones', 'rand')
p.addParameter('initialProbe','circ')   % choose initial probe ('obes', 'rand')
p.addParameter('intensityConstraint','standard') % ('standard', 'sigmoid')
p.addParameter('makeGIF', false) % probe orthogonalization frequency
p.addParameter('npsm',1)                % number of probe state mixtures
p.addParameter('nosm',1)                % number of object state mixtures
p.addParameter('numIterations',1)       % number of iterations
p.addParameter('objectPlot','complex')  % 'complex', 'abs', 'angle', piAngle
p.addParameter('objectUpdateStart',1)   % iteration when object update starts (makes algorithm stable if good initial guess is known)
p.addParameter('orthogonalizationFrequency', 10) % probe orthogonalization frequency
p.addParameter('positionOrder','random')% position order for sequential solvers ('sequential' or 'random')
p.addParameter('probePowerCorrectionSwitch', true)  % switch to indicate if probe power correction is applied
p.addParameter('saveMemory', true)      % if true, then the algorithm is run with low memory consumption (e.g. detector error is not saved)
p.addParameter('singleSwitch', true)    % if true, data is converted to single to save memory and processing time
p.addParameter('smoothenessSwitch', false)  % if true, apply smoothness constraint to object
p.addParameter('objectTVregSwitch', false) % TV regularization
p.addParameter('gpuSwitch', false)      % if true, then algorithm runs on gpu
p.addParameter('zPIEgradientStepSize', 100) % gradient step size for axial position correction (typical range [1, 100])

p.parse( varargin{:} )

% engine
obj.params.engine = p.Results.engine;

% ePIE parameters
obj.params.betaObject = p.Results.betaObject;
obj.params.betaProbe = p.Results.betaProbe;

% algorithmic parameters
obj.params.absorbingProbeBoundary = p.Results.absorbingProbeBoundary;
obj.params.absObjectSwitch = p.Results.absObjectSwitch;
obj.params.positionOrder = p.Results.positionOrder;
obj.params.intensityConstraint = p.Results.intensityConstraint;
obj.params.probePowerCorrectionSwitch = p.Results.probePowerCorrectionSwitch;
obj.params.smoothenessSwitch = p.Results.smoothenessSwitch;
obj.params.comStabilizationSwitch = p.Results.comStabilizationSwitch;
obj.params.gpuSwitch = p.Results.gpuSwitch;
obj.params.saveMemory = p.Results.saveMemory;
obj.params.orthogonalizationFrequency = p.Results.orthogonalizationFrequency;
obj.params.batchSize = p.Results.batchSize;
obj.params.zPIEgradientStepSize = p.Results.zPIEgradientStepSize;
obj.params.objectTVregSwitch = p.Results.objectTVregSwitch;
% general parameters
obj.params.npsm = p.Results.npsm;
obj.params.nosm = p.Results.nosm;

% detector
obj.params.deleteFrames = p.Results.deleteFrames;
obj.params.FourierMaskSwitch = p.Results.FourierMaskSwitch;

% reconstruction
obj.params.numIterations = p.Results.numIterations;
obj.params.figureUpdateFrequency = p.Results.figureUpdateFrequency;
obj.params.objectPlot = p.Results.objectPlot;

% processing
obj.params.singleSwitch = p.Results.singleSwitch;
obj.params.fftshiftSwitch = p.Results.fftshiftSwitch;
obj.params.makeGIF = p.Results.makeGIF;
obj.monitor.fontSize = p.Results.fontSize;

% initialize detector error matrices
if obj.params.saveMemory
    obj.params.detectorError = 0;
else
    obj.params.detectorError = zeros(obj.Nd, obj.Nd, obj.numFrames);
end

if ~isempty(obj.ptychogram)
    obj.params.energyAtPos = squeeze(sum(sum( abs(obj.ptychogram), 1 ), 2));
else
    obj.params.energyAtPos = squeeze(sum(sum( abs(obj.params.ptychogramDownsampled), 1 ), 2));
end

% detector (delete Frames if necessary)
if ~isempty(obj.params.deleteFrames)
    
    if  ~isempty(obj.ptychogram)
        obj.ptychogram(:,:,obj.params.deleteFrames) = [];
    else
        obj.params.ptychogramDownsampled(:,:,obj.params.deleteFrames) = [];
    end
    
    if ~obj.params.saveMemory
        obj.params.detectorError(:,:,obj.params.deleteFrames) = [];
    end
    obj.params.energyAtPos(obj.params.deleteFrames) = [];
    
    if ~isempty(obj.positions)
        obj.positions(obj.params.deleteFrames,:) = [];
    end
    
    if ~isempty(obj.positions0)
        obj.positions0(obj.params.deleteFrames,:) = [];
    end
    
    try
        obj.params.encoder(obj.params.deleteFrames,:) = [];
    catch
        disp('encoder not set in raw data')
    end
    
    if ~isempty(obj.ptychogram)
        obj.numFrames = size(obj.ptychogram, 3);
    else
        obj.numFrames = size(obj.params.ptychogramDownsampled, 3);
    end
    
end

%% save memory

if obj.params.saveMemory
    obj.Xo = [];
    obj.Yo = [];
end

% initial guesses

% initial object
switch p.Results.initialObject
    case 'ones'
        obj.params.initialObject = ones(obj.No, obj.No, obj.params.nosm) + ...
            0.001 * rand(obj.No, obj.No, obj.params.nosm);
        
    case 'rand'
        obj.params.initialObject = exp(1i * pi * (rand(obj.No, obj.No, obj.params.nosm)-1/2));
        obj.params.initialObject = convolve2( obj.params.initialObject, gaussian2D(5,3), 'same' );
end

% initialProbe
switch p.Results.initialProbe
    
    case {'circ', 'Circ', 'CIRC'}
        obj.params.initialProbe = ones(obj.Np, obj.Np, obj.params.npsm);
        
        pupil = circ(obj.Xp, obj.Yp, obj.entrancePupilDiameter);
        obj.params.initialProbe = bsxfun(@times, obj.params.initialProbe, pupil);
        
        % force linear independency of initial guess if there are multiple (state mixture) modes
        obj.params.initialProbe = obj.params.initialProbe .* (0.99 + 0.01 * rand(obj.Np, obj.Np, obj.params.npsm));
        
    case {'rand', 'Rand', 'RAND'}
        obj.params.initialProbe = exp(1i*2*pi*(rand(obj.Np, obj.Np, obj.params.npsm)));
        
        % for CPM set initial guess to entrance pupil (e.g. pinhole diameter)
        pupil = circ(obj.Xp, obj.Yp, obj.entrancePupilDiameter);
        obj.params.initialProbe = bsxfun(@times, obj.params.initialProbe, pupil);
        
    case {'gaussian', 'Gaussian', 'GAUSSIAN'}
        stdGauss = obj.entrancePupilDiameter / 2 / (2.355);
        gaussianBeam = exp(-(obj.Xp.^2+obj.Yp.^2)/(2*stdGauss^2));
        obj.params.initialProbe = [];
        for k = 1:(obj.params.npsm)
            obj.params.initialProbe = cat(3, obj.params.initialProbe, gaussianBeam);
        end
        
        % randomize to force linear independence
        if obj.params.npsm > 1
            obj.params.initialProbe = obj.params.initialProbe .* (0.99 + 0.01 * rand(obj.Np, obj.Np, obj.params.npsm));
        end
        
    case 'initialGuess'
        disp('use initial guess contained in obj.params.initialProbe')
    otherwise
        error('obj.params.initialProbe incorrectly set.')
        
end

% quadratic phase Term
if strcmp(obj.propagator.type, 'Fraunhofer')
    obj.propagator.quadraticPhase = single(exp(1i * pi/(obj.wavelength*obj.zo) * (obj.Xp.^2 + obj.Yp.^2)));
    obj.params.initialProbe = bsxfun(@times, obj.params.initialProbe, obj.propagator.quadraticPhase);
elseif strcmp(obj.propagator.type, 'Fresnel')
    obj.propagator.quadraticPhase = single(exp(1i * pi/(obj.wavelength*obj.zo) * (obj.Xp.^2 + obj.Yp.^2)));
    obj.params.initialProbe = bsxfun(@times, obj.params.initialProbe, obj.propagator.quadraticPhase);
end


%% initial OTF (only for TRAM reconstruction)

% absorbing probe boundary
if or(not(obj.params.saveMemory), obj.params.absorbingProbeBoundary)
    % filter probe with super-gaussian window function
    c = 4/5; m  = 4;
    obj.params.probeWindow = exp(-( (obj.Xp.^2)/(2*(c*obj.Np*obj.dxp/2.355)^2 ) ).^m) .* ...
        exp(-( (obj.Yp.^2)/(2*(c*obj.Np*obj.dxp/2.355)^2 ) ).^m);
end

% normalize probe to energy in ptychogram
if ~isempty(obj.ptychogram)
    obj.params.probePowerCorrection = sqrt( max( sum(sum(obj.ptychogram, 1), 2)) );
else
    obj.params.probePowerCorrection = sqrt( max( sum(sum(obj.params.ptychogramDownsampled, 1), 2)) );
end
obj.params.initialProbe = obj.params.initialProbe / sqrt(sum(abs(obj.params.initialProbe(:)).^2)) * obj.params.probePowerCorrection;

% initialize class properties

% set estimated positions equal to initial positions
if isempty(obj.positions)
    obj.positions = obj.positions0;
end

if isempty(obj.numFrames)
    if ~isempty(obj.ptychogram)
        obj.numFrames = size(obj.ptychogram, 3);
    else
        obj.numFrames = size(obj.params.ptychogramDownsampled, 3);
    end
end

if isempty(obj.object)
    obj.object = obj.params.initialObject;
    if obj.params.saveMemory
        obj.params.initialObject = [ ];
    end
end

if isempty(obj.probe)
    obj.probe = obj.params.initialProbe;
    if obj.params.saveMemory
        obj.params.initialProbe = [ ];
    end
end

%% initialize properties

% set "if-exist" quantities (quantities that should not be overwritten)
if ~isfield(obj.params, 'error')
    % initialize error
    obj.params.error = [];
end

if ~isfield(obj.params ,'objectROI')
    % set object ROI (for monitoring)
    %     obj.params.objectROI = [min(obj.positions0(:,1))+obj.Np/2, max(obj.positions0(:,1))+obj.Np/2,...
    %         min(obj.positions0(:,2))+obj.Np/2, max(obj.positions0(:,2))+obj.Np/2];
    obj.params.objectROI = [min(obj.positions0(:,1)), max(obj.positions0(:,1))+obj.Np,...
        min(obj.positions0(:,2)), max(obj.positions0(:,2))+obj.Np];
    
end

if ~isfield(obj.params, 'probeROI')
    % set probe ROI
    n = 1;
    r = obj.entrancePupilDiameter/obj.dxp;
    obj.params.probeROI = ...
        [max(1,obj.Np/2 - n*r), min(obj.Np,obj.Np/2 + n*r),...
        max(1,obj.Np/2 - n*r), min(obj.Np,obj.Np/2 + n*r)];
    obj.params.probeROI = round(obj.params.probeROI);
    
end

if ~isfield(obj.params, 'errorAtPos')
    obj.params.errorAtPos = zeros(obj.numFrames, 1);
end

if ~isfield(obj.params, 'runtime')
    obj.params.runtime = tic;
end

if ~isfield(obj.params,'singleSwitch')
    % initialize single switch to false if not already set elsewhere
    obj.params.singleSwitch = false;
end

if ~isfield(obj.params,'verboseLevel')
    % initialize verbose level as 'low' if not already set elsewhere
    obj.params.verboseLevel = 'low';
end

if obj.params.singleSwitch
    obj.convert2single;
end

if ~isfield(obj.params,'fftshiftFlag')
    % normally the raw data should not be shifted, so the flag is set to 0
    obj.params.fftshiftFlag = 0;
end

if ~isfield(obj.params, 'gpuFlag')
    obj.params.gpuFlag = 0;
end

if ~isfield(obj.params, 'CPSCswitch')
    obj.params.CPSCswitch = false;
end

%% monitor

if ~isfield(obj.monitor, 'DIFFcmap')
    obj.monitor.DIFFcmap = setColormap;
end

%% calculate detector NA and expected depth of field

obj.params.NAd = obj.Ld/(2*obj.zo);
obj.params.DoF = obj.wavelength(1) / obj.params.NAd^2;

return
