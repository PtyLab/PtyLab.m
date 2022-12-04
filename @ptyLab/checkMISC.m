function checkMISC(obj)
% checkMISC:        checks miscellaneous quantities specific certain engines

% shuffle random number generator
rng('shuffle')

if obj.params.backgroundModeSwitch
   obj.params.background = ones(obj.Np, obj.Np, 'like', obj.probe); 
end

% preallocate intensity scaling vector
if strcmp(obj.params.intensityConstraint, 'fluctuation')
   obj.params.intensityScaling = ones(1, obj.numFrames, 'like', obj.params.error);
end

% check if there is data on gpu that shouldn't be there
if and( isfield(obj.params, 'diffractionData'), ~strcmp( obj.params.engine, 'kPIE' ) )
    obj.params.diffractionData = [ ];
end

% check if both probe power correction and modulus enforced probe are
% switched on. Since this can cause a contradiction, it raises an erro
if and(obj.params.modulusEnforcedProbeSwitch, obj.params.probePowerCorrectionSwitch)
    error('Modulus enforced probe and probe power correction can not simultaneaously be switched on!')
end

if strcmp( obj.propagator.type, 'ASP' )
    
    [~, obj.propagator.transferFunction]  = aspw(obj.probe(:,:,1), obj.zo, obj.lambda, obj.Lp);
    obj.propagator.transferFunction = single(obj.propagator.transferFunction);
    
    if obj.params.fftshiftSwitch
       error('ASP propagator works only with obj.params.fftshiftSwitch = false!') 
    end
end

if strcmp( obj.propagator.type, 'scaledASP' )
    
    [~, obj.ODpropagator.Q1, obj.ODpropagator.Q2]  = scaledASP(obj.probe(:,:,1), obj.zo, obj.lambda, obj.dxo, obj.dxd);
    
    if obj.params.fftshiftSwitch
       error('scaledASP propagator works only with obj.params.fftshiftSwitch = false!') 
    end
end

if ~obj.params.fftshiftSwitch
    warning('fftshiftSwitch set to false; this may lead to reduced performance')
end