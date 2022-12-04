function checkDataset(obj)
% checkDataset:     examines potential issues in the data

% check positions0
if size(obj.positions0,2) > 2
    error('Property "positions0" is set in the wrong way. See class definition for help.')
end

% check if ptychogram (or at least downsampled version is available)
if isempty(obj.ptychogram)
    if or( ~isfield(obj.params,'ptychogramDownsampled'), isempty(obj.params.ptychogramDownsampled) )
        error('ptychogram (or downsampled version) needs to be supplied')
    elseif isfield(obj.params,'ptychogramDownsampled')
        obj.params.CPSCswitch = true;
    end
end


disp('===============================')
disp('checkDataset: no problems found')
disp('===============================')
return