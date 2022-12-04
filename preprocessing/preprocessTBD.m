% preprocessing
 
clear
restoredefaultpath
tic

% add toolbox folders
switch getenv('COMPUTERNAME')
    case 'PHASESPACE'
        toolboxFolder = 'F:\Dropbox\Codes\ptyLab'; addpath(genpath(toolboxFolder)); 
        obj = ptyLab; obj.export.exportPath = 'F:\ptyLabExport';
        obj.export.filePath = '\\sun\fkoenderink\shared-projects\stefanfemius-shared-projects\AntennaPtycho\Data_Ptymore\201124_noAntennaTest\20201124_144147_USAF2\AVT camera (GT2300)';
    otherwise
        toolboxFolder = '/home/user/Dropbox/Codes/ptyLab'; addpath(genpath(toolboxFolder));
        obj = ptyLab; obj.export.exportPath = '/home/user/Documents/ptyLabExport';
        obj.export.filePath = '--';
        
end
cd(toolboxFolder);

% change this depending on computer used for preprocessing

addpath(genpath(obj.export.filePath))

obj.export.fileName = 'recent';  % this can be defined and determines output file name
camera = 'GT2300';          % 'Hamamatsu', 'GX'
% wavelength
obj.wavelength = 632.8e-9;   % blue laser diode
% binning
obj.binningFactor = 2;
% padding for superresolution
padFactor = 1;
% set magnification if any objective lens is used
magnification = 1;
% object detector distance (initial guess)
obj.zo = 31.3e-3;
% dark / readout offset

% set verbose level ('high' or 'low' depending on how many details should be shown)
obj.params.verboseLevel = 'low';
% set detection geometry
% A: camera to closer side of stage (allows to bring camera close in transmission)
% B: camera to further side of stage (doesn't allow to bring close in transmission)
% ...other way around in reflection
% C: objective + tube lens in transmission
measurementMode = 'A';

%%

switch camera
    case 'Hamamatsu'
        M = 2^11; % Hammamatsu
        N = 2^11; % Hammamatsu
        obj.dxd = 6.5e-6 * obj.binningFactor / magnification;  % effective detector pixel size is magnified by binning
        backgroundOffset = 30;      % globally subtracted from raw data (diffraction intensities), play with this value
    case 'GX'
        N = 1456; % GX
        M = 1456; % GX
        obj.dxd = 4.54e-6 * obj.binningFactor / magnification;  % effective detector pixel size is magnified by binning
        backgroundOffset = 30;       % globally subtracted from raw data (diffraction intensities), play with this value (good value for GX)
    case 'GT2300'
        M = 1752;
        N = 1752;
        obj.dxd = 5.5e-6 * obj.binningFactor / magnification;
        backgroundOffset = 30;
end

%% number of frames is calculated automatically 

% number of Frames
temp = dir(obj.export.filePath);
obj.numFrames = length(temp) - 5; % this assumes the folder only containes 

%% read background

dark = single(imread('background.tif'));

%% read empty beam (if available)

try 
    emptyBeam = single(imread('emptyBeam.tif'));
    emptyBeamBool = true;
catch
    emptyBeamBool = false;
end
    

%% binning

obj.ptychogram = zeros(N/obj.binningFactor * padFactor, ...
                       N/obj.binningFactor * padFactor, obj.numFrames, 'single');

disp('==== execute binning ====')
for k = 1:obj.numFrames
    % get file name
    if k<10
        fileNameTemp = ['0000',num2str(k-1),'.tif'];
    elseif and( 9 < k-1, k-1 < 100)
        fileNameTemp = ['000',num2str(k-1),'.tif'];
    elseif and( 99 < k-1, k-1 < 1000)
        fileNameTemp = ['00',num2str(k-1),'.tif'];
    elseif and( 999 < k-1, k-1 < 10000)
        fileNameTemp = ['0',num2str(k-1),'.tif'];
    end
    disp(['...reading frame ',fileNameTemp])
    % read images and correct for background
    temp = posit(single(imread(fileNameTemp)) - dark - backgroundOffset * ones(size(dark),'single') );
    
    % crop
%     temp = temp(M/2-N/2+1:M/2+N/2, M/2-N/2+1:M/2+N/2, :);
    temp = temp(end/2-M/2+1:end/2+M/2, end/2-N/2+1:end/2+N/2, :);
    
     % flipping and binning
     switch measurementMode
         case 'A'
             temp = flipud( imresize( temp , 1/obj.binningFactor, 'nearest' ) );
         case 'B'
             temp = rot90( imresize( temp, 1/obj.binningFactor, 'nearest' ), 2 );
         case 'C'
             temp = rot90(flipud( imresize( temp , 1/obj.binningFactor, 'nearest' ) ), 2);
     end
    
     %zero padding
     [obj.ptychogram(:,:,k), W] = diff_mask( temp, N/obj.binningFactor*padFactor * [1,1]);

end
disp(['binning time:', num2str(toc)])

if emptyBeamBool
    switch measurementMode
        case 'A'
            obj.params.emptyBeam = diff_mask( flipud( imresize( emptyBeam , 1/obj.binningFactor, 'nearest' ) ),  N/obj.binningFactor*padFactor * [1,1]);
        case 'B'
            obj.params.emptyBeam = diff_mask( rot90( imresize( emptyBeam , 1/obj.binningFactor, 'nearest' ), 2 ),  N/obj.binningFactor*padFactor * [1,1]);
        case 'C'
            error('case C todo')
    end
end

% note: sometimes, the data has to be flipped, because the camera
% manufacturer outputs the data as it would appear either in or against
% beam propagation direction; in some setups the camera is rotated/flipped,
% so always check this

%% set experimental specifications

% detector coordinates
obj.Nd = size(obj.ptychogram, 1);       % number of detector pixels
obj.Ld = obj.Nd * obj.dxd;              % effective size of detector
obj.xd = (-obj.Nd/2:obj.Nd/2-1)*obj.dxd;% 1D coordinates in detector plane
[obj.Xd, obj.Yd] = meshgrid(obj.xd);    % 2D coordinates in detector plane

% object coordinates
obj.dxo = obj.wavelength * obj.zo / obj.Ld;   % Fraunhofer/Fresnel
% obj.dxo = obj.dxd;                      % asp
% obj.dxo = 400e-9 * obj.zo / obj.Ld;       % scaled asp, choose freely (be careful not to depart too much from Fraunhofer condition)
obj.No = 2^12;
obj.Lo = obj.No * obj.dxo;
obj.xo = (-obj.No/2:obj.No/2-1) * obj.dxo;
[obj.Xo, obj.Yo] = meshgrid(obj.xo);

% probe coordinates
obj.dxp = obj.dxo;
obj.Np = obj.Nd;
obj.Lp = obj.Np * obj.dxp;
obj.xp = (-obj.Np/2:obj.Np/2-1)*obj.dxp;
[obj.Xp, obj.Yp] = meshgrid(obj.xp);

%% Fourier mask
% variant 1
% obj.params.W = W;
% obj.params.W(~W) = 0.001;
% n = 20;
% obj.params.W(1:n,:) = 1;
% obj.params.W(end-n:end,:) = 1;
% obj.params.W(:,1:n) = 1;
% obj.params.W(:,end-n:end) = 1;

obj.params.W = circ(obj.Xd, obj.Yd, obj.Nd/2 * obj.dxd) > 0;
obj.params.W(~obj.params.W) = 0.001;

% variant 2
% obj.params.W = 1-exp(-(obj.Xd.^2+obj.Yd.^2)/(1.25*obj.Nd*obj.dxd/2.355).^2);
% obj.params.W(obj.params.W>1) = 1;
% obj.params.W(W == 1) = 1;
% 
% figure(99)
% imagesc(obj.params.W)
% axis image
% colorbar

%% set entrance pupil diameter 

obj.entrancePupilDiameter = 1000e-6;

%% get positions

% get file name (this assumes there is only one text file in the raw data folder)
positionsFileName = dir([obj.export.filePath,'\*.txt']);


% 2) take raw data positions
T = readtable(positionsFileName.name,'delimiter', ' ');
% T = readtable(positionsFileName.name);
try
    positions0 = [T.y, T.x];
catch
    positions0 = [T.Var1, T.Var2];
end
    
positions0(:,1) = positions0(:,1) - positions0(1,1);
positions0(:,2) = positions0(:,2) - positions0(1,2);

% center positions
% positions0 = bsxfun(@minus, positions0, positions0(1,:));
% convert positions to micrometer
% positions0 = positions0 * 1e-6 / magnification;
positions0 = positions0 * 1e-6;

obj.params.encoder = positions0;
% convert positions into pixels
positions0 = round( positions0 / obj.dxo );
% center within object grid
positions0 = positions0 + obj.No/2 - obj.Np/2;
% take only the frames needed (if numFrames smaller than the number of positions in the file)
obj.positions0 = positions0(1:obj.numFrames,:); 

% show positions
figure(3)
plot(positions0(:,2), positions0(:,1),'o-')
axis equal,axis image
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

exportBool = true;
if exportBool
%     obj.exportID = [obj.fileName,'_bin',num2str(obj.binningFactor)];
    obj.export.exportID = obj.export.fileName;
    obj.exportObj
end
toc

 