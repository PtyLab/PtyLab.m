% preprocessing
 
clear
restoredefaultpath
tic


% add relevant folders, specified in config.m
s = pwd;
switch s(end-5:end)
    case 'ptyLab'
        addpath(genpath('config'))
    otherwise
        addpath(genpath('../config'))     
end
addpath(genpath('../config'));
[toolboxFolder, dataFolder] = config;
        
addpath(genpath( toolboxFolder ));
addpath(genpath( dataFolder ));

cd(toolboxFolder);

% create ptyLab object
obj = ptyLab('operationMode', 'CPM');
obj.export.exportPath = dataFolder;

% change this depending on computer used for preprocessing
obj.export.filePath = 'F:\Dropbox\CurrentProjects\ptychoWorkshop\data\lenspaper4\AVT camera (GX1920)';
addpath(genpath(obj.filePath))

obj.fileName = 'recent';  % this can be defined and determines output file name
camera = 'GX';          % 'Hamamatsu', 'GX'
% wavelength
obj.wavelength = 450e-9;   % NIR
% binning
obj.binningFactor = 2;
% padding for superresolution
padFactor = 1;
% set magnification if any objective lens is used
magnification = 1;
% object detector distance (initial guess)
obj.zo = 19.31e-3;
% dark / readout offset

% set verbose level ('high' or 'low' depending on how many details should be shown)
obj.params.verboseLevel = 'low';
% set detection geometry
% A: camera to closer side of stage (allows to bring camera close in transmission)
% B: camera to further side of stage (doesn't allow to bring close in transmission)
% ...other way around in reflection
% C: objective + tube lens in transmission
measurementMode = 'A';

%% correct for miscalibration of AOTF using look-up table

obj.lambda = AOTFwavelengthLUT(obj.lambda);

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
end

%% number of frames is calculated automatically 

% number of Frames
temp = dir(obj.filePath);
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
    fileNameTemp = strcat(num2str(k-1, '%05d'), '.tif');
    
    disp(['...reading frame ',fileNameTemp])
    % read images and correct for background
    temp = posit(single(imread(fileNameTemp)) - dark - backgroundOffset * ones(size(dark),'single') );
    
    % crop
    temp = temp(M/2-N/2+1:M/2+N/2, M/2-N/2+1:M/2+N/2, :);
    
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
%% set power spectral density

obj.params.PSD = mean( obj.ptychogram, 3);

%% set experimental specifications

% detector coordinates
obj.Nd = size(obj.ptychogram, 1);       % number of detector pixels
obj.Ld = obj.Nd * obj.dxd;              % effective size of detector
obj.xd = (-obj.Nd/2:obj.Nd/2-1)*obj.dxd;% 1D coordinates in detector plane
[obj.Xd, obj.Yd] = meshgrid(obj.xd);    % 2D coordinates in detector plane

% object coordinates
obj.dxo = obj.lambda * obj.zo / obj.Ld;   % Fraunhofer/Fresnel
% obj.dxo = obj.dxd;                      % asp
% obj.dxo = 400e-9 * obj.zo / obj.Ld;     % scaled asp, choose freely 
% note: (be careful not to depart too much from Fraunhofer condition)
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

%% high NA (curvature) correction
% % obj.zo = 1e-2;
% Xn = obj.Xd / (obj.Ld/2);
% Yn = obj.Yd / (obj.Ld/2);
% NA = obj.Ld/(2*obj.zo);
% rp = sqrt(Xn.^2 + Yn.^2);
% Xq = Xn ./ sqrt(1 + (rp * NA).^2);
% Yq = Yn ./ sqrt(1 + (rp * NA).^2);
% 
% if NA > 0.1
%     disp('curvature correction...')
%     correction1 = 1 + (obj.Xd.^2 + obj.Yd.^2)/obj.zo;
%     
%     for k = 1:obj.numFrames
%         obj.ptychogram(:,:,k) = obj.ptychogram(:,:,k) .* correction1;
%         obj.ptychogram(:,:,k) = posit(interp2( Xn, Yn, obj.ptychogram(:,:,k),Xq, Yq,'cubic',0));
% %         figure(99)
% %         %     imagesc(log10(obj.ptychogram(:,:,k)+1))
% %   
% %         axis image, colormap jet
% %         drawnow
% %         pause(0.5)
% %         title(num2str(k))
%     end
%     
%     W = interp2( Xn, Yn, W,Xq, Yq,'cubic',0) > 0.9;
% end
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
positionsFileName = dir([obj.filePath,'\*.txt']);

% 2) take raw data positions
T = readtable(positionsFileName.name);
positions0 = [T.Var1, T.Var2];
positions0(:,1) = positions0(:,1) - positions0(1,1);
positions0(:,2) = positions0(:,2) - positions0(1,2);

% center positions
% positions0 = bsxfun(@minus, positions0, positions0(1,:));
% convert positions to micrometer
% positions0 = positions0 * 1e-6 / magnification;
positions0 = positions0 * 1e-6;

obj.params.encoder = positions0;
% convert positions into pixels
positions0 = positions0 / obj.dxo;
% center within object grid
positions0 = positions0 + obj.No/2 - obj.Np/2;
% then round to nearest integer
positions0 = round(positions0);
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
    obj.exportID = obj.fileName;%[obj.fileName,'_bin',num2str(obj.binningFactor)];
    obj.exportObj
end
toc

 