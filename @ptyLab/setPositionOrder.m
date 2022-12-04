function indices = setPositionOrder(obj)

switch obj.params.positionOrder
    case 'sequential'
        indices = 1:obj.numFrames;
    case 'random'
        if isempty(obj.params.error)
            indices = 1:obj.numFrames;
        else
            if length(obj.params.error) < 2
                indices = 1:obj.numFrames;
            else
                indices = randperm(obj.numFrames);
            end
        end
    case 'FP' % brightfield first
        rows = obj.positions(:,1) - mean( obj.positions(:,1) );
        cols = obj.positions(:,2) - mean( obj.positions(:,2) );
        dist = sqrt(rows.^2 + cols.^2);
        [~, indices] = sort(dist, 'ascend');
    otherwise
        error('position order not properly set')
end