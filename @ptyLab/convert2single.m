function obj = convert2single( obj )
% function to convert to single precision
%
% use this function to accelerate data processing
% note: be aware of reduced numerical precision of single as compared to
% float variables

display('==== convert to single ====')

if strcmp(obj.params.verboseLevel,'high')
    
    display('object size before converting to single precision')
    objSize = obj.getSize('showPieChart', true,'figureNumber', 98);
    title(gca,{'memory distribution before conversion to single',...
        ['total memory: ', num2str(round(objSize/1e6)),' MB']} )
    
end

% convert properties to single precision
if isprop(obj,'diffData')
    obj.diffData = single(obj.diffData);
end

if isprop(obj,'ptychogram')
    obj.ptychogram = single(obj.ptychogram);
end

if isfield(obj.params,'PSD')
    obj.params.PSD = single(obj.params.PSD);
end

if isprop(obj,'object')
    obj.object = single(obj.object);
end

if isprop(obj,'probe')
    obj.probe = single(obj.probe);
end

if isprop(obj,'Xo')
    obj.Xo = single(obj.Xo);
end

if isprop(obj,'Yo')
    obj.Yo = single(obj.Yo);
end

if isprop(obj,'Xd')
    obj.Xd = single(obj.Xd);
end

if isprop(obj,'Yd')
    obj.Yd = single(obj.Yd);
end

if isprop(obj,'Xp')
    obj.Xp = single(obj.Xp);
end

if isprop(obj,'Yp')
    obj.Yp = single(obj.Yp);
end

% convert parameters to single precision
if isfield(obj.params,'initialObject')
    obj.params.initialObject = single(obj.params.initialObject);
end

if isfield(obj.params, 'detectorError')
   obj.params.detectorError = single(obj.params.detectorError); 
end

if isfield(obj.params, 'probeWindow')
   obj.params.probeWindow = single(obj.params.probeWindow); 
end

% convert ODpropagator to single precision

if isfield(obj.propagator, 'QitmIn')
   obj.propagator.QitmIn = single(obj.ODpropagator.QitmIn); 
end

if isfield(obj.propagator, 'QitmOut')
   obj.propagator.QitmOut = single(obj.ODpropagator.QitmOut); 
end

if isfield(obj.propagator, 'QIn')
   obj.propagator.QIn = single(obj.ODpropagator.QIn); 
end

if isfield(obj.propagator, 'QOut')
   obj.propagator.QOut = single(obj.ODpropagator.QOut); 
end

if strcmp(obj.params.verboseLevel,'high')
    display('object size after converting to single precision')
    
    %% display size and show pie chart
    objSize = obj.getSize('showPieChart', true,'figureNumber',99);
    title(gca,{'memory distribution after conversion to single',...
        ['total memory: ', num2str(round(objSize/1e6)),' MB']} )
end

end