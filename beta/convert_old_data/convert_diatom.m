addpath(genpath('C:\Users\Lars Loetgering\Dropbox\Codes\fracPty'))
clear
load('C:\Users\Lars Loetgering\Documents\fracPtyExport\MPI_190706021_recLarge_MPIstar.mat')

%% write hdf5 file

try
    delete *.h5
end

file_name = 'C:\Users\Lars Loetgering\Documents\fracPtyRaw\ptyLab_data\ptyLab_helical_beam.h5';
% file_name = 'C:\Users\Lars Loetgering\Documents\fracPtyRaw\ptyLab_data\ptyLab_OAM_spokes.h5';
disp('write h5 file data')
tic
% write ptychogram
hdf5write(file_name, ...
         '/ptychogram', obj.ptychogram(:,:,idx),... 
         '/encoder', obj.positions0(idx,:)*obj.dxo,...
         '/dxd', obj.dxd,...
         '/zo', obj.zo,...
         '/wavelength', 2.48e-9);

ptychogram = h5read(file_name,'/ptychogram');
params.encoder = h5read(file_name,'/encoder')';
dxd = h5read(file_name,'/dxd');
zo = h5read(file_name,'/zo');
wavelength = h5read(file_name,'/wavelength');
numFrames = size(obj.ptychogram,3);
toc


