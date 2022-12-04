% MAXYMUS preprocessing

close all
clear
tic

% add relevant folders
toolboxFolder = 'C:\Users\Lars Loetgering\Dropbox\Codes\PtyLab'; 
addpath(genpath(toolboxFolder)); obj = ptyLab;
file_name = 'C:\Users\Lars Loetgering\Dropbox\fracPtyRaw';
cd(toolboxFolder);

obj.export.fileName = 'MPI_180623067';
obj.export.exportPath = 'C:\Users\Lars Loetgering\Documents\fracPtyExport';
% % get export path
% obj.export.getExportPath 

% create data link (data is not read yet)
fid = H5F.open(['C:\Users\Lars Loetgering\Dropbox\fracPtyRaw\MPI_180623067.cxi']);

% get intensities
M = 0;
% N = 2^8 + M;
N = 2^8 - 2^4;
obj.binningFactor = 1;

%% read data

dset_id = H5D.open(fid,'/entry_1/data_1/data');
intensity_patterns = H5D.read(dset_id);  % reads e.g. N x N x nop diffraction patterns (intensities)

switch obj.export.fileName
    
    case 'MPI_180623059'
        centerRow = 133;
        centerCol = 132;
    case 'MPI_180623064'
        centerRow = 132;
        centerCol = 132;
    case 'MPI_180623066' % helical beam paper, data set 1
        centerRow = 132;
        centerCol = 132;
    case 'MPI_180623067'
        centerRow = 132;
        centerCol = 132;
    case 'MPI_190118006'
        centerRow = 130;
        centerCol = 130;
    case 'MPI_190118011'
        centerRow = 130;
        centerCol = 130;
    case 'MPI_190118012'
        centerRow = 130;
        centerCol = 130;
    case 'MPI_190118026'
        centerRow = 129;
        centerCol = 126;
    case 'MPI_190118030'
        centerRow = 129;
        centerCol = 126;
    case 'MPI_190119008'
        centerRow = 132;
        centerCol = 132;
    case 'MPI_190119013'
        centerRow = 132;
        centerCol = 132;
    case 'MPI_190119011'
        centerRow = 133;
        centerCol = 132;
    case 'MPI_190119019'
        centerRow = 133;
        centerCol = 132;
    case 'MPI_190120019'
        centerRow = 132;
        centerCol = 131;
    case 'MPI_190120021'
        centerRow = 132;
        centerCol = 131;
    otherwise
        error('center not explicitly specified yet for this filename!')
end

intensity_patterns = intensity_patterns(centerRow-N/2+1:centerRow+N/2,centerCol-N/2+1:centerCol+N/2,:);
obj.numFrames = size(intensity_patterns, 3);
%% enforce positivity on data

intensity_patterns = posit(intensity_patterns);

%%% note: 
%%% try to not set all this data to zero but write a decent background
%%% correction algorithm
% intensity_patterns = intensity_patterns - min(intensity_patterns(:));

%% binning

if obj.binningFactor > 1
    % execute binning
    obj.ptychogram = zeros(N/obj.binningFactor, N/obj.binningFactor, obj.numFrames, 'single');
    
    for k = 1:obj.numFrames
        obj.ptychogram(:,:,k) = imresize( intensity_patterns(:,:,k) , 1 / obj.binningFactor, 'nearest' ) ;
    end
    
else
    % do not execute binning
    obj.ptychogram = intensity_patterns;
end

%% get power spectral density

obj.params.PSD = mean( obj.ptychogram, 3 );

%%
% convert to single precision to save memory
obj.ptychogram = single(obj.ptychogram); 
% clear intensity_patterns

% get wavelength
energy_id = H5D.open(fid,'/entry_1/instrument_1/source_1/energy');
energy = H5D.read(energy_id)/(1.602e-19);
obj.wavelength = 1240 / energy * 1e-9;

% get detector pixel size
x_pixel_size = H5D.open(fid,'/entry_1/instrument_1/detector_1/x_pixel_size');
obj.dxd = H5D.read(x_pixel_size) * obj.binningFactor; % detector pixel size (pixel size quadratic)
obj.Nd = size(obj.ptychogram, 1);
obj.Ld = obj.Nd * obj.dxd;
obj.xd = (-obj.Nd/2:obj.Nd/2-1)*obj.dxd;
[obj.Xd, obj.Yd] = meshgrid(obj.xd);

% get object detector distance
switch obj.export.fileName
    case 'MPI_180623066' % helical beam paper, data set 1
        obj.zo = 169e-3;
    case 'MPI_190118006'
        obj.zo = 21e-2;
    case 'MPI_190118011'
        obj.zo = 21e-2;
    case 'MPI_190118012'
        obj.zo = 21e-2;
    case 'MPI_190118026'
        obj.zo = 21e-2;
    case 'MPI_190118030'
        obj.zo = 21e-2;
    case 'MPI_190119011'
        obj.zo = 22e-2;
    case 'MPI_190119008'
        obj.zo = 21e-2;
    case 'MPI_190119013'
        obj.zo = 20e-2;
    case 'MPI_190120019'
        obj.zo = 20.5e-2;
    case 'MPI_190120021'
        obj.zo = 20.5e-2;
    
    otherwise
        z_id = H5D.open(fid,'/entry_1/instrument_1/detector_1/distance');
        obj.zo = H5D.read(z_id);
end
obj.dxo = obj.wavelength * obj.zo / obj.Ld;
obj.No = 2^12;
obj.Lo = obj.No * obj.dxo;
obj.xo = (-obj.No/2:obj.No/2-1)*obj.dxo;
[obj.Xo, obj.Yo] = meshgrid(obj.xo);

% probe coordinates
obj.dxp = obj.dxo;
obj.Np = obj.Nd;
obj.Lp = obj.Np * obj.dxp;
obj.xp = (-obj.Np/2:obj.Np/2-1)*obj.dxp;
[obj.Xp, obj.Yp] = meshgrid(obj.xp);

% get positions
translation_id = H5D.open(fid,'/entry_1/sample_1/geometry_1/translation');
translations = H5D.read(translation_id);
obj.positions0 = [translations(1,:)' translations(2,:)'];   % 1) [Y, X] (use row column order)
obj.positions0 = bsxfun(@minus, obj.positions0, sum(obj.positions0, 1)/obj.numFrames );   % 2) center positions
% obj.positions0 = bsxfun(@minus, obj.positions0, obj.positions0(1,:) );   % 2) center positions
obj.positions0 = obj.positions0 * 1;        % this line is a reminder that the MAXYMUS positions are in units of meters and do not need to be converted
obj.positions0 = round(obj.positions0 / obj.dxo);   % convert to pixel units
% obj.positions0 = obj.positions0 + obj.No/2 - obj.Np/2;  % center within object grid
obj.positions = obj.positions0;     % copy 'positions0' into position estimate 'positions'

% set object-detector propagator
obj.propagator.type = 'Fraunhofer';

% probe size was 1µm in most experiments
obj.entrancePupilDiameter = 1.5e-6;

%% define Fourier mask


obj.params.W = 1-circ(obj.Xd + obj.dxd, obj.Yd + obj.dxd, obj.dxd * 4);

% obj.params.W(62+M/2:66+M/2,:) = 0; % apply only if centerRow = 133; centerCol = 131; !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% obj.params.W(126,124:131) = 0;
% obj.params.W(127,125) = 0;
% 
% % line in first quadrant
% obj.params.W(1:end/2,114+M/2) = 0;
% 
% % defect pixels insecond quadrant
% obj.params.W(60+M/2:63+M/2,1+M/2:4+M/2) = 0;
% 
% % strange signal in third quadrant
% obj.params.W(end/2:end,12+M/2:13+M/2) = 0;
% 
% % blob (4th quadrant)
% obj.params.W(97+M/2,126+M/2:127+M/2) = 0;
% obj.params.W(98+M/2,125+M/2:128+M/2) = 0;  
% obj.params.W(99+M/2,126+M/2:127+M/2) = 0;
% 
% % circular region
% obj.params.W = obj.params.W .* (1-circ(obj.Xd + obj.dxd, obj.Yd + obj.dxd, obj.dxd * 5));
% obj.params.W = obj.params.W > 0;
% obj.params.W0 = obj.params.W;

figure(10)
hsvplot( log10( obj.params.PSD + 1 ) .* exp(1i * obj.params.W) )
colormap jet
axis image off

%% export data

exportBool = true;
if exportBool
    obj.export.exportID = [obj.export.fileName,'_bin',num2str(obj.binningFactor)];
    obj.exportObj
end

%% close hdf5 file

H5D.close(dset_id);
H5F.close(fid);

%% clear IDs
clear dset_id energy_id energy z_id x_pixel_size translation_id fid
toc

%% get only every second row from raw data

idx = [ ];

% % way 1: small FOV
for k = 1:80
    if mod(k,2) == 0
    else
        idx = [idx, ((k-1)*40+(1:40))];
    end
end


% way 2: checkerboard moderate FOV
% for k = 1:80
%     if mod(k,2) == 0
%         idx = [idx, ((k-1)*80+(1:2:80))];
%     else
%         idx = [idx, ((k-1)*80+(2:2:80))];
%     end
% end

file_name = 'C:\Users\Lars Loetgering\Documents\fracPtyRaw\ptyLab_data\ptyLab_helical_beam.h5';
disp('write h5 file data')
tic
% write ptychogram
hdf5write(file_name, ...
         '/ptychogram', obj.ptychogram(:,:,idx),... 
         '/encoder', obj.positions0(idx,:)*obj.dxo,...
         '/dxd', obj.dxd,...
         '/zo', obj.zo,...
         '/wavelength', 2.48e-9);

