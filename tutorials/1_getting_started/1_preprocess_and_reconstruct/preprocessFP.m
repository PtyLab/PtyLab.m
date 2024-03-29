% preprocess FP data set
 
clear
restoredefaultpath
tic

% add toolbox folders
directory_info = dir();
toolboxFolder = erase( directory_info(1).folder, 'tutorials\1_getting_started\1_preprocess_and_reconstruct');

cd(toolboxFolder);
addpath(genpath(toolboxFolder))

obj = ptyLab;
% change this depending on computer used for preprocessing
obj.export.exportPath = 'C:\Users\Lars Loetgering\Documents\ptyLabExport';
obj.export.filePath = 'C:\Users\Lars Loetgering\Documents\fracPtyRaw\ptyLab_data\LungCarcinomaFPM.hdf5';
% obj.export.filePath = 'C:\Users\Lars Loetgering\Documents\fracPtyRaw\ptyLab_data\USAFTargetFPM.hdf5';
obj.export.fileName = 'recent'; % this determines output file name
obj.params.verboseLevel = 'low';
obj.operationMode = 'FPM';

%% read hdf5 file

tic
disp('read h5 file data')
obj.ptychogram = single(h5read(obj.export.filePath,'/ptychogram'));
obj.params.encoder = h5read(obj.export.filePath,'/encoder')';
obj.dxd = h5read(obj.export.filePath,'/dxd');
obj.wavelength = h5read(obj.export.filePath,'/wavelength');
obj.params.magnification = h5read(obj.export.filePath,'/magnification');
obj.params.zled = h5read(obj.export.filePath,'/zled');
obj.params.NA = h5read(obj.export.filePath,'/NA');
obj.numFrames = size(obj.ptychogram,3);
toc

%% correct orientation of encoder and ptychogram

obj.ptychogram = fliplr(rot90( obj.ptychogram, 1));
obj.params.encoder = -[obj.params.encoder(:,1), obj.params.encoder(:,2)];

%% set experimental specifications

% detector coordinates
obj.Nd = size(obj.ptychogram, 1);       % number of detector pixels
obj.Ld = obj.Nd * obj.dxd;              % effective size of detector
obj.xd = (-obj.Nd/2:obj.Nd/2-1)*obj.dxd;% 1D coordinates in detector plane
[obj.Xd, obj.Yd] = meshgrid(obj.xd);    % 2D coordinates in detector plane

% object spectrum coordinates (or lens plane)
obj.dxo = 1 / (obj.Ld / obj.params.magnification);   % FPM
obj.No = 2^9;          % this may have to be adjusted, depending on data set
obj.Lo = obj.No * obj.dxo;
obj.xo = (-obj.No/2:obj.No/2-1) * obj.dxo;
[obj.Xo, obj.Yo] = meshgrid(obj.xo);

% pupil coordinates
obj.dxp = obj.dxo;
obj.Np = obj.Nd;
obj.Lp = obj.Np * obj.dxp;
obj.xp = (-obj.Np/2:obj.Np/2-1)*obj.dxp;
[obj.Xp, obj.Yp] = meshgrid(obj.xp);

%% convert k-space positions (as given by obj.params.encoder) into integer-valued positions
% convert positions into pixels
prefactor = obj.Ld / obj.params.magnification / obj.wavelength;
obj.positions0 = round( prefactor * obj.params.encoder ./ ...
    sqrt(obj.params.encoder(:,1).^2 + obj.params.encoder(:,2).^2 + obj.params.zled^2) );
% center within object grid
% obj.positions0  = obj.positions0 - floor(mean(obj.positions0, 1))  + obj.No/2 - round(obj.Np/2);
% obj.positions0  = obj.positions0 - floor(mean(obj.positions0, 1))  + obj.No/2;
obj.positions0  = obj.positions0 + obj.No/2 - round(obj.Np/2);

% take only the frames needed (if numFrames smaller than the number of positions in the file)
obj.positions0  = obj.positions0(1:obj.numFrames,:);

% show positions
figure(3)
plot(obj.positions0(:,2), obj.positions0(:,1),'o-')
axis equal, axis image
set(gcf, 'color', 'w')

%% set propagator

obj.propagator.type = 'Fraunhofer';

%% save memory by converting to single precision 
% (saves roughly a factor of 1.5 to 2 in memory and reconstruction speed)

obj.params.singleSwitch = true;
if obj.params.singleSwitch
   obj.convert2single; 
end

%% export data

obj.export.format = 'mat';
exportBool = true;
if exportBool
    obj.export.exportID = obj.export.fileName;
    obj.exportObj
end
toc

 