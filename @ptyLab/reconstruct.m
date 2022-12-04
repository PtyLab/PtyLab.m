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
    case 'OPRP'
        % orthogonal probe relaxation (OPR)
        obj = OPRP(obj);
    case 'pcPIE'
        % position correctiion PIE
        obj = pcPIE(obj);   
    case 'e3PIE'
        % multiSlice PIE
        obj = e3PIE(obj);
    case 'zPIE' 
        % axial position correction
        obj = zPIE(obj);
    case 'pSD'
        % parallel steepest descent
        obj = pSD(obj);
    case 'PIE'
        % original PIE algorithm
        obj = PIE(obj);
    case 'adam_FP'
        % accelerated solver for FP
        obj = adam_FP(obj);
    case 'adam'
        % accelerated solver for CP
        obj = adam(obj);
    case 'qNewton'
        % quasi-Newton algorithm
        obj = qNewton(obj);
    case 'mqNewton'
        % momentum-accelerated quasi Newton
        obj = mqNewton(obj);
    otherwise
        error('engine not specified')
        
end

end