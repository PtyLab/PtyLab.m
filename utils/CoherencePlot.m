function [W,mu] = CoherencePlot(probe, r, direction)
% example: 

row = r(1);
col = r(2);

% regularization
aleph = 1e-3;

% preallocation
W = zeros( size(probe, 1), size(probe, 2) );

switch direction
    case 'x'
    
        for k = 1:size(probe, 3)
            W = W + probe(row,:,k)'*probe(row,:,k);
            wdiag = sqrt(diag(W));
            denom = wdiag * wdiag';
            mu = W./(denom + aleph);
        end
        
    case 'y'
        for k = 1:size(probe, 3)
            W = W + probe(:,col,k)*probe(:,col,k)';
            wdiag = sqrt(diag(W));
            denom = wdiag * wdiag';
            mu = W./(denom + aleph);
        end
end

% renormalization to compensate regularization
mu = mu/max(abs(mu(:)));

end