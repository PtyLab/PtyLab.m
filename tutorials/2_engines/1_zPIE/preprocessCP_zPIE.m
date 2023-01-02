% preprocessing

clear
restoredefaultpath
tic

set(0,'DefaultFigureWindowStyle','normal')
set(0, 'DefaultAxesFontSize', 18);
set(0, 'DefaultFigureColor', 'w');
set(0,'defaultAxesFontName', 'serif')

% add toolbox folders
directory_info = dir();
toolboxFolder = erase( directory_info(1).folder, 'tutorials\2_engines\1_zPIE');

%%
cd(toolboxFolder);
addpath(genpath(toolboxFolder))

obj = ptyLab;
% change this depending on computer used for preprocessing
obj.export.exportPath = 'C:\Users\Lars Loetgering\Documents\ptyLabExport';
obj.export.filePath = 'C:\Users\Lars Loetgering\Documents\fracPtyRaw\ptyLab_data\ptyLab_zPIE.h5';
obj.export.fileName = 'recent'; % this determines output file name
obj.params.verboseLevel = 'low';

%% read hdf5 file

tic
disp('read h5 file data')
obj.ptychogram = single(h5read(obj.export.filePath,'/ptychogram'));
obj.params.encoder = h5read(obj.export.filePath,'/encoder');
obj.dxd = h5read(obj.export.filePath,'/dxd');
obj.zo = h5read(obj.export.filePath,'/zo');
obj.wavelength = h5read(obj.export.filePath,'/wavelength');
obj.numFrames = size(obj.ptychogram,3);
toc

%% set z-distance incorrectly 
% this will be corrected during the reconstruction using the zPIE
% calibration engine

obj.zo = 33e-3;

%% rotate ptychogram and encoder to monitor data in upright position

obj.ptychogram = fliplr(rot90(obj.ptychogram, 1));
obj.params.encoder = [-obj.params.encoder(:,2), -obj.params.encoder(:,1)];

%% sampling

% detector coordinates
obj.Nd = size(obj.ptychogram, 1);
obj.Ld = obj.Nd * obj.dxd;
obj.xd = (-obj.Nd/2:obj.Nd/2-1) * obj.dxd;
[obj.Xd, obj.Yd] = meshgrid(obj.xd);    % 2D coordinates in detector plane

% object coordinates
obj.dxo = obj.wavelength * obj.zo / obj.Ld;   % Fraunhofer/Fresnel
obj.No = 2^11;
obj.Lo = obj.No * obj.dxo;
obj.xo = (-obj.No/2:obj.No/2-1) * obj.dxo;
[obj.Xo, obj.Yo] = meshgrid(obj.xo);

% probe coordinates
obj.dxp = obj.dxo;
obj.Np = obj.Nd;
obj.Lp = obj.Np * obj.dxp;
obj.xp = (-obj.Np/2:obj.Np/2-1)*obj.dxp;
[obj.Xp, obj.Yp] = meshgrid(obj.xp);

%% convert positions

% convert positions into pixels
obj.positions0 = round( (obj.params.encoder) / obj.dxo );
% center within object grid
obj.positions0  = obj.positions0 - floor(mean(obj.positions0, 1))  + obj.No/2 - round(obj.Np/2);
% take only the frames needed (if numFrames smaller than the number of positions in the file)
obj.positions0  = obj.positions0(1:obj.numFrames,:);

% show positions
figure(3)
plot(obj.positions0 (:,2), obj.positions0(:,1),'o-')
axis equal,axis image
set(gcf, 'color', 'w')

%% set entrance pupil diameter

obj.entrancePupilDiameter = obj.Lp/3; % rough beam size estimate
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
    obj.export.exportID = [obj.export.fileName];
    obj.exportObj
end

%% show raw data
obj.showPtychogram
toc

