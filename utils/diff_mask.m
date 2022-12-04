function [psi_d, S]=diff_mask(diff_data,varargin)
% [psi_d S]=diff_mask(diff_data,varargin)
% alternatively: [psi_d S]=diff_mask(diff_data,[M N])
%
% diff_data: gray scale diffraction pattern
%
% diff_mask() returns diffraction pattern embedded within (2^M x 2^N)-image, where N
% is the smallest integer such that psi_d is larger than diff_data
%
% S is the a binary mask to access the diffraction data
%
% author: Lars Loetgering (lars.loetgering@fulbrightmail.org)
% last change: 09/26/14

[M, N] = size(diff_data);

switch nargin
    case 1
        Mp = 2^(ceil(log2(M)));
        Np = 2^(ceil(log2(N)));
    case 2
        A = varargin{1};
        Mp = A(1); 
        Np=A(2);
    otherwise;
        error('wrong number of input arguments')
end

S = zeros(Mp,Np,'like', diff_data);
psi_d = zeros(Mp,Np,'like', diff_data);
dM = Mp-M;
dN = Np-N;

if mod(dM,2)~=0
    error('To produce even number of output pixels both dimensions of input should be even.')
elseif mod(dN,2)~=0
    error('To produce even number of output pixels both dimensions of input should be even.')
end

S(dM/2+1:dM/2+M,dN/2+1:dN/2+N) = 1;
psi_d(dM/2+1:dM/2+M,dN/2+1:dN/2+N) = diff_data;

end