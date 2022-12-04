function checkMemory(obj)
% checkMemory:      keeps the memory consumption low by removing fields
% that are not needed for the reconstruction engine specified

if obj.params.saveMemory
    
    if and( isfield(obj.params, 'probeStack'), ~strcmp(obj.params.engine, 'OPRP') )
        obj.params = rmfield(obj.params, 'probeStack');
        obj.params = rmfield(obj.params, 'OPRPmodes');
        obj.params = rmfield(obj.params, 'OPRPpurity');
        obj.params = rmfield(obj.params, 'normalizedOPRPEigenvalues');
        obj.params = rmfield(obj.params, 'OPRPU');
        obj.params = rmfield(obj.params, 'OPRPS');
        obj.params = rmfield(obj.params, 'OPRPV');
    end
    
    if ~strcmp(obj.params.engine, 'mPIE')
        if isfield(obj.params, 'objectBuffer')
            obj.params = rmfield(obj.params, 'objectBuffer');
        end
        
        if isfield(obj.params, 'objectMomentum')
            obj.params = rmfield(obj.params, 'objectMomentum');
        end
        
        if isfield(obj.params, 'probeBuffer')
            obj.params = rmfield(obj.params, 'probeBuffer');
        end
        
        if isfield(obj.params, 'probeMomentum')
            obj.params = rmfield(obj.params, 'probeMomentum');
        end
    end
end

