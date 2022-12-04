function [probe, normalizedEigenvalues, V] = orthogonalizeModes(probe, varargin)
% [probe, normalizedEigenvalues] = orthogonalizeModes(probe, numModes)
% imposes orthogonality through singular value decomposition

if nargin > 1
    numModes = varargin{1};
else
    numModes = size(probe,3);
end

p = reshape(probe, [size(probe,1)*size(probe,2), size(probe,3)]);

% method of snapshots
[V, D] = eig( p' * p );
S = real(sqrt(D));
U = (p * V ) / (S + eps * eye(size(D)));

% sort in descending order
S = diag(S);
[Ssorted,idx] = sort(S, 'descend');
S = diag(Ssorted');
U = U(:,idx);
V = V(:,idx);
probe = reshape(U*S, size(probe,1), size(probe,2), numModes);
normalizedEigenvalues = gather( diag(S.^2)/trace(S.^2) );

% truncated SVD method
% [U, S, V] = fsvd(p, numModes, 1, 1);
% probe = reshape(U*S, size(probe,1), size(probe,2), numModes);
% normalizedEigenvalues = gather( diag(S.^2)/trace(S.^2) );


% %% Gram-Schmidt
% 
% n = size(p,1);
% U = zeros(n,numModes,'like',p);
% U(:,1) = p(:,1)/sqrt(p(:,1)'*p(:,1));
% 
% for i = 2:numModes
%     U(:,i) = p(:,i);
%     for j = 1:i-1
%         U(:,i) = U(:,i) - ( U(:,j)'*U(:,i) ) / ( U(:,j)'*U(:,j) ) * U(:,j);
%     end
%     normalizedEigenvalues = (U(:,i)'*U(:,i));
% %     U(:,i) = U(:,i)/sqrt(U(:,i)'*U(:,i));
% end
% 
% probe = reshape(U, size(probe,1), size(probe,2), numModes);
end