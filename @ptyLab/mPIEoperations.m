function mPIEoperations( obj, varargin )
% mPIEoperations:   Several extra operations for mPIE are contained in this function

p = inputParser;
p.addParameter('mode', 'init')
p.parse( varargin{:} );

switch p.Results.mode
    case 'init'
        % initialize momentum and time (if necessary)
%         if and(~isfield(obj.params, 'objectMomentum'), strcmp(obj.params.engine, 'mPIE'))
        if ~isfield(obj.params, 'objectMomentum')
            obj.params.objectMomentum = zeros(size(obj.object),'like', obj.object);
        elseif and(~isfield(obj.params, 'objectSlicesMomentum'), strcmp(obj.params.engine, 'm3PIE'))
            obj.params.objectSlicesMomentum = zeros(size(obj.params.objectSlices),'like', obj.params.objectSlices);
        else
            % TODO
            % allow for dynamic object mode structure
%             if obj.params.nosm < size( obj.params.objectMomentum, 3 )
%                 error('dynamic object mode structure not yet implemented for mPIE')
%             elseif obj.params.nosm > size( obj.params.objectMomentum, 3 )
%                 error('dynamic object mode structure not yet implemented for mPIE')
%             end
        end
        
        if ~isfield(obj.params, 'probeMomentum')
            obj.params.probeMomentum = zeros(size(obj.probe),'like', obj.probe);
        else
            % allow dynamic mode structure (if npsm/nosm changes)
            if obj.params.npsm < size( obj.params.probeMomentum, 3 )
                obj.params.probeMomentum = obj.params.probeMomentum(:,:,1:obj.params.npsm);
            elseif obj.params.npsm > size( obj.params.probeMomentum, 3 )
                numModesAdded =  obj.params.npsm - size( obj.params.probeMomentum, 3 );
                obj.params.probeMomentum =cat(3, obj.params.probeMomentum, zeros(obj.Np, obj.Np, numModesAdded, 'like', obj.probe) );
            end
            
        end
        
    case 'orthogonalize'
        obj.orthogonalize
        
        if obj.params.npsm > 1
            % orthogonalize probe Buffer
            p = reshape(obj.params.probeBuffer, [obj.Np^2, obj.params.npsm]);
            obj.params.probeBuffer = reshape(p * obj.params.MSPVprobe, size(obj.probe));
            % orthogonalize probe momentum
            p = reshape(obj.params.probeMomentum, [obj.Np^2, obj.params.npsm]);
            obj.params.probeMomentum = reshape(p * obj.params.MSPVprobe, size(obj.probe));
        end
        
    case 'clearMemory'
        % this case is only executed at the end of an mPIE reconstruction
        obj.params = rmfield(obj.params, 'probeBuffer');
        
    otherwise
        error('error in mPIEhcecks: mode not specified')
end

end