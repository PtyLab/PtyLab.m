function exportObj(obj)

if isunix
    fileName = [obj.export.exportPath,'/',obj.export.exportID,'.mat'];
else
    fileName = [obj.export.exportPath,'\',obj.export.exportID,'.mat'];
end

save(fileName,'obj','-v7.3')
disp([fileName,' exported'])
        
end
