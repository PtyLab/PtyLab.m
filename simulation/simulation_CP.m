% configuration
set(0, 'DefaultAxesFontSize', 18);
set(0, 'DefaultFigureColor', 'w');
set(0,'defaultAxesFontName', 'serif')
clear
close all
restoredefaultpath

% set the toolbox folder in the next line 
% (i.e. the folder where you saved PtyLab)
toolboxFolder = 'C:\Users\Lars Loetgering\Dropbox\Codes\PtyLab';
addpath(genpath(toolboxFolder))
% set the data folder in the next line 
% (i.e. the folder where you plan to export the data 
% generated in this simulation)
dataFolder  = 'C:\Users\Lars Loetgering\Documents\ptyLabExport';
addpath(dataFolder)

%% create PtyLab object

obj = ptyLab('operationMode', 'CPM'); % generate object (of class PtyLab)
% note: type open 'PtyLab' to see the properties of this class
% - in the class definition you find required fills which are currently
% empty 
% - the purpose of this code is to fill in those properties and simulate 
% a CPM data set. To this end, all properties need to be defined. 

exportBool = true;              % if exportBool = true, data is saved in dataFolder
% get export path
obj.export.exportPath = dataFolder;    % this is where simulated data is saved

%% set physical properties

obj.wavelength = 632.8e-9; % wavelength
obj.zo = 5e-2;             % object detector distance

% detector coodinates (in practice, sampling requirements are set by 
% the detector! So start from here)
obj.Nd = 2^7;                       % number of samples in detector plane
obj.dxd = 2^11 / obj.Nd * 4.5e-6;   % detector pixel size
obj.Ld = obj.Nd * obj.dxd;          % detector size
obj.xd = (-obj.Nd/2:obj.Nd/2-1)*obj.dxd;% 1D coordinates (detector)
[obj.Xd, obj.Yd] = meshgrid(obj.xd);    % 2D coordinates (detector)

% entrance pupil/ probe coordinates
obj.dxp = obj.wavelength * obj.zo / obj.Ld; % probe sampling step size
obj.Np = obj.Nd;                % number of samples in probe field of view
obj.Lp = obj.Np * obj.dxp;      % field of view in pinhole plane
obj.xp = (-obj.Np/2:obj.Np/2-1)*obj.dxp;% 1D coordinates (probe)
[obj.Xp, obj.Yp] = meshgrid(obj.xp);% 2D coordinates (probe)
obj.zp = 1e-2;                  % pinhole-object distance

% object coordinates
obj.dxo = obj.dxp;              % object pixel size
obj.No = 2^10 + 2^9;            % number of samples in object plane
obj.Lo = obj.No * obj.dxo;      % object field of view
obj.xo = (-obj.No/2:obj.No/2-1)*obj.dxo;% 1D coordinates (object)
[obj.Xo, obj.Yo] = meshgrid(obj.xo);% 2D coordinates (object)

%% generate illumination
% note: simulate focused beam 
% goal: 1:1 image iris through (low-NA) lens with focal length f onto an object

f = 2.5e-3;  % focal length of lens, creating a focused probe
pinhole = circ(obj.Xp, obj.Yp, obj.Lp/2); % type "help circ" to see the definition of the circ function
% note: the below uncommented lines may be used to vary the illumination profile
% pinhole = pinhole .* single(imread('cameraman.tif')); 
pinhole = pinhole .* exp(1i * atan2(obj.Yp,obj.Xp)); % OAM beam
pinhole = normconv2( pinhole, gaussian2D(5,1) );
% note: the convolution of the circular aperture makes the edges of the lens smooth 

% propagate to lens
probe = aspw(pinhole, 2*f, obj.wavelength, obj.Lp);  

% multiply with quadratic phase and aperture
aperture = circ(obj.Xp, obj.Yp, 3*obj.Lp/4);
aperture = normconv2( aperture, gaussian2D(5,3) );
% apply lens
probe = probe .* exp(-1i * 2*pi/obj.wavelength * (obj.Xp.^2 + obj.Yp.^2)/(2*f)) .* aperture;

% probe = aspw(probe, 2*f, obj.wavelength, obj.Lp); 
probe = aspw(probe, 2*f, obj.wavelength, obj.Lp); 
% note: aspw = angular spectrum of plane waves (propagates field right behind lens by distance 2f)

% show probe (complex amplitude)
figure(1); 
subplot(1,2,1)
hsvplot(obj.xp, probe);  
% note: color = phase, brightness = intensity [more precisely, brightness = sqrt(intensity) )
colormap gray; axis image, title('complex probe')
% show probe (intensity)
subplot(1,2,2)
imagesc(obj.xp,obj.xp, abs(probe).^2);
axis image, colormap gray
title('probe intensity')
%% generate object

d = 1e-3;   % the smaller this parameter the larger the spatial frequencies in the simulated object
b = 33;     % topological charge (feel free to play with this number)
[theta, rho] = cart2pol(obj.Xo, obj.Yo);
t = ( 1 + sign( sin( b * theta + 2*pi * (rho/d).^2) ) )/2;
% phaseFun = exp(1i * atan2(obj.Yo, obj.Xo));
phaseFun = 1;
% phaseFun = exp(1i*( 1 * theta + 2*pi * (rho/d).^2));
t = t .* circ(obj.Xo,obj.Yo, obj.Lo) .* (1 - circ(obj.Xo,obj.Yo,200*obj.dxo)) .* phaseFun + ...
    circ(obj.Xo,obj.Yo,130*obj.dxo);
object = normconv2( t, gaussian2D(5,3) ); % smooth edges

% load object
figure(2); hsvxplot(object,'pixelSize',obj.dxo, 'colorbar', 'test')
axis image, title('object')

%% generate positions

% scan parameters
numPoints = 200;            % number of points
radius = 100;               % radius of final scan grid [in micrometers]
p = 1;                      % "clumping parameter"* (see note below)
% expected beam size, required to calculate overlap (expect Gaussian-like beam, derive from second moment)
beamSize = sqrt(sum2((obj.Xp.^2 + obj.Yp.^2) .* abs(probe).^2) / sum2(abs(probe).^2)) * 2.355; 
% * note:
% * p = 1 is standard Fermat 
% * p > 1 yields more points towards the center of grid   

% generate non-uniform Fermat grid
[R, C] = GenerateNonUniformFermat(numPoints, 'radius', radius, 'power', p);

% 2 opt
[R,C] = twoOpt(R,C);

% prevent negative indices bz centering spiral coordinates on object
R = R + size(object, 1)/2 - obj.Np/2 + 50;
C = C + size(object, 2)/2 - obj.Np/2 + 50;
R = round(R); C = round(C);

% get number of positions
obj.numFrames = length(R); 
disp(['generate positions (',num2str(obj.numFrames),')'])

% show scan grid
figure(99)
plot(R, C, 'ko-','MarkerFaceColor','k')
axis square
title('scan grid')
axis([min(C), max(C), min(R), max(R)])

%% generate ptychogram

obj.ptychogram = zeros(obj.Nd, obj.Nd, obj.numFrames);
cmap = setColormap;
% specify propagator between sample and detector (Fraunhofer, Fresnel, ASP, scaledASP)
obj.propagator.type = 'Fresnel'; 
obj.propagator.quadraticPhase = ...
    exp(1i * 2*pi/obj.wavelength * (obj.Xp.^2 + obj.Yp.^2)/(2*obj.zo));
obj.params.fftshiftSwitch = false;
for loop = 1:obj.numFrames
    
    % get object patch
    objectPatch = object( R(loop) : R(loop) + obj.Np-1, C(loop) : C(loop)+obj.Np-1 );
    
    % multiply each probe mode with object patch
    obj.params.esw = probe .* objectPatch;
    
    % generate diffraction data (complex amplitude conversion below)
    obj.object2detector;
    
    % save data in ptychogram 
    I = sum( abs(obj.params.ESW).^2, 3 );
    obj.ptychogram(:, :, loop) = I;
    
    % inspect diffraction data
    figure(10)
    imagesc(obj.xd, obj.xd, log10(I + 1e-3)); axis image; colormap(cmap); colorbar; 
    title('diffraction data')
    drawnow
end

%%  calculate noise

% simulate Poisson noise
bitDepth = 12;
maxNumCountsPerDiff = 2^bitDepth;

% normalize data (ptychogram)
I = I/max(obj.ptychogram(:)) * maxNumCountsPerDiff;
obj.ptychogram = obj.ptychogram/max(obj.ptychogram(:)) * maxNumCountsPerDiff;

% simulate Poisson noise
obj.ptychogram(:) = poissrnd(obj.ptychogram(:)); 

% compare noiseless data to noisy data
figure(12)
imagesc(sqrt([I, obj.ptychogram(:,:,loop)])) 
axis image off; colormap(cmap) 
title(['left: noiseless, right: noisy (',num2str(bitDepth),' bit)'])

%% set data 

% note: the entrance pupil diameter determines the probe/pupil initial estimate diameter 
obj.entrancePupilDiameter = beamSize;
obj.positions0 = [R,C];
obj.params.encoder = (obj.positions0 - obj.No/2 + obj.Np/2) * obj.dxo;

%% data inspection, check sampling requirements

obj.checkDataset

%% create slider through raw data

obj.showPtychogram;

%% export data

obj.export.format = 'mat';
if exportBool
    obj.export.exportID = 'recent';
    obj.exportObj
end
