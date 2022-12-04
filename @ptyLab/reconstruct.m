function obj = reconstruct(obj)
% start measuring time (toc is used in engines)
tic

% check quantities (mode structure, fftshifts etc.)
obj.checkMemory
obj.checkModes
obj.checkFFT
obj.checkHandles
obj.checkMISC
obj.checkGPU
disp('start reconstruction')

switch obj.params.engine
    
    % steepest descent
    case 'SD'
        obj = SD(obj); % sequential SD
    % PIE-type engines
    case 'ePIE' 
        % standard ePIE
        obj = ePIE(obj);
    case 'lsqPIE'
        % adaptive step size PIE 
        obj = lsqPIE(obj);
    case 'mPIE' 
        % momentum PIE
        obj = mPIE(obj);
    case 'pcPIE'
        % position correctiion PIE
        obj = pcPIE(obj); 
    case 'OPRP'
        % orthogonal probe relaxation (OPR)
        obj = OPRP(obj);
    case 'e3PIE'
        % multiSlice PIE
        obj = e3PIE(obj);
    case 'm3PIE'
        % multiSlicePIE with momentum
        obj = m3PIE(obj);
    case 'kPIE' % using k means clustering
        obj = kPIE(obj);
    case 'zPIE' % axial position correction
        obj = zPIE(obj);
    case 'pSD'
        obj = pSD(obj);
    case 'PIE'
        obj = PIE(obj);
    case 'adam_FP'
        obj = adam_FP(obj);
    case 'adam'
        obj = adam(obj);
    case 'qNewton'
        obj = qNewton(obj);
    case 'mqNewton'
        obj = mqNewton(obj);
    otherwise
        error('engine not specified')
end

end