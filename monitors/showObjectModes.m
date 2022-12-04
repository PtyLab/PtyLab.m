function showObjectModes(obj)

hsvmodeplot( obj.object(obj.params.objectROI(1):obj.params.objectROI(2), ...
    obj.params.objectROI(3):obj.params.objectROI(4),: ) )

end
    