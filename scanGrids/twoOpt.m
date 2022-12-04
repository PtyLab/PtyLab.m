function [R,C] = twoOpt(R, C, varargin)
% heuristic optimizer for traveling salesman problem 
% R: row coordinates
% C: column coordinates
% > each should be a column vector
% source: https://github.com/larsloetgering/2opt

p = inputParser;
p.addParameter('closed_path', true) % constant multiplier
p.parse(varargin{:})

if size(R,1) < size(R,2)
    R = R'; C = C';
end
disp('--')
disp('optimize travel path')
% nearest-neighbor preconditioning
tic
route_precon = precondition_route(R,C);
t1 = toc;
disp(['nearest neighbor preconditioning: ', num2str(t1), ' sec'])

% 2opt algorithm
tic
route2Opt = two_opt(route_precon(:,1),route_precon(:,2));
t2 = toc;
disp(['2opt: ', num2str(t2), ' sec'])
disp('--')

R = route2Opt(:,1);
C = route2Opt(:,2);

if p.Results.closed_path
    % makes sure the first scan point is the same as the last scan point
    R = [R; R(1)]; C = [C; C(1)];
end

end

%% PRECONDITIONING
function route = precondition_route(R,C)
numPos = size(R,1);
idx = 2:numPos-1;
idx_nearest = [1;zeros(numPos-2, 1);numPos];
for k=2:numPos-1
    distances = (R(idx_nearest(k-1))-R(idx)).^2 + (C(idx_nearest(k-1))-C(idx)).^2;
    next_idx = find(distances == min(distances));
    idx_nearest(k) = idx(next_idx(1));
    idx = setdiff(idx, idx(next_idx(1)));
end
R = R(idx_nearest);
C = C(idx_nearest);
route = [R,C];
end
%% 2OPT
function route = two_opt(R,C)

route = [R,C];
best_route = route;
cost_best = cost(best_route);
progress_made = true;
while progress_made
    progress_made = false;
    for k = 1:length(route)-1
        for l=k+1:length(route)
            % reverse potential crossing
            new_route = [route(1:k,:); ...
                        flipud(route(k+1:l-1,:));...
                        route(l:end,:)];
            % update best route (if lowest cost achieved)
            cost_new = cost(new_route);
            if cost_new < cost_best
                best_route = new_route;
                cost_best = cost_new;
                progress_made = true;
            end
            
        end
    end
    route = best_route;
end
end
%% COST
function r = cost(route)
r = sum(diff(route(:,1)).^2 + diff(route(:,2)).^2);
end