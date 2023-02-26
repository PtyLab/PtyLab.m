% note: the purpose of this simulation script is to illustrate the
% conversion of conventional ptychography (CP) data to Fourier ptychography
% (FP) data

% configuration

clear
close all
set(0, 'DefaultAxesFontSize', 18);
set(0, 'DefaultFigureColor', 'w');
set(0,'defaultAxesFontName', 'serif')

%% add toolbox and data folder

directory_info = dir();
toolboxFolder = erase( directory_info(1).folder, 'tutorials\5_conversion_CP_FP');
% adjust the export path to your local file system
exportPath  = 'C:\Users\Lars Loetgering\Documents\ptyLabExport';

cd(toolboxFolder) 
addpath(genpath( toolboxFolder ));  % toolbox contains all required functions to run code

%% create ptyLab object

obj = ptyLab('operationMode', 'CPM'); % generate object (of class ptyLab)
% note: type open 'ptyLab' to see the properties of this class
% - in the class definition you find required fields which are currently
% empty 
% - the purpose of this code is to fill in those properties and simulate 
% a CPM data set. To this end, all properties need to be defined.
% - subsequently, we will convert what originally used to be a CP data set
% into a FP data set, thus illustrating the conversion from CP to FP

exportBool = true;  % if exportBool = true, data is saved in dataFolder
% get export path
obj.export.exportPath = exportPath;% this is where simulated data is saved
clear exportPath % not needed anymore, so clean up

%% set physical properties

obj.wavelength = 632.8e-9;      % wavelength
obj.zo = 5e-2;                  % object detector distance

% detector coodinates (typically sampling requirements are set by the detector! So start from here)
obj.dxd = 2^11 / 2^6 * 4.5e-6;  % detector pixel size (this is equivalent to 8-fold binning)
obj.Nd = 2^6;                   % number of samples in detector plane
obj.Ld = obj.Nd * obj.dxd;      % detector size
obj.xd = (-obj.Nd/2:obj.Nd/2-1)*obj.dxd;% 1D coordinates (detector)
[obj.Xd, obj.Yd] = meshgrid(obj.xd);    % 2D coordinates (detector)

% entrance pupil/ probe coordinates
% here we use the sampling for conventional Fresnel/Fraunhofer ptychography
obj.dxp = obj.wavelength * obj.zo / obj.Ld; % probe sampling step size
obj.Np = obj.Nd;                % number of samples in probe field of view
obj.Lp = obj.Np * obj.dxp;      % field of view in pinhole plane
obj.xp = (-obj.Np/2:obj.Np/2-1)*obj.dxp;% 1D coordinates (probe)
[obj.Xp, obj.Yp] = meshgrid(obj.xp);% 2D coordinates (probe)
obj.zp = 1e-2;                  % pinhole-object distance

%% generate illumination
% note: simulate focused beam 
% goal: 1:1 image pinhole through (low-NA) lens with focal length f onto an object

f = 2e-3;  % focal length of lens, creating a focused probe
pinhole = circ(obj.Xp, obj.Yp, obj.Lp/2); % type "help circ" to see the definition of the circ function
pinhole = normconv2( pinhole, gaussian2D(5,1) ); % smooth edges
% note: the convolution of the circular aperture makes the edges of the lens smooth 

% propagate to lens
probe = aspw(pinhole, 2*f, obj.wavelength, obj.Lp);  

% multiply with quadratic phase (=lens) and aperture
aperture = circ(obj.Xp, obj.Yp, 3*obj.Lp/4);
aperture = normconv2( aperture, gaussian2D(5,3) ); % smooth edges
probe = probe .* exp(-1i * 2*pi/obj.wavelength * ...
                     (obj.Xp.^2 + obj.Yp.^2)/(2*f)) .* aperture;

probe = aspw(probe, 2*f, obj.wavelength, obj.Lp);
% note: aspw = angular spectrum of plane waves 
% (propagates field right behind lens by distance 2f)

% show focused probe
figure(1); 
subplot(1,2,1)
hsvplot(probe);  
% note: color = phase, brightness = intensity [more precisely, brightness = sqrt(intensity) )
colormap gray; axis image, title('complex probe')

subplot(1,2,2)
imagesc(obj.xp,obj.xp, abs(probe).^2);
axis image, colormap gray
title('probe intensity')

%% generate positions
S = 32; % number of scan points per dimension
dS = 4;
beamSize = sqrt(sum2((obj.Xp.^2 + obj.Yp.^2) .* abs(probe).^2) / sum2(abs(probe).^2)) * 2.355; 
beamSizePix = beamSize / obj.dxp;

[R, C] = generate_raster_grid(S, dS);
R = round(R); C = round(C);

% get number of positions
obj.numFrames = length(R); 
disp(['generate positions (',num2str(obj.numFrames),')'])

%% generate object

% object coordinates
obj.dxo = obj.dxp;              % object pixel size
obj.No = 512;
obj.Lo = obj.No * obj.dxo;      % object field of view
obj.xo = (-obj.No/2:obj.No/2-1)*obj.dxo;% 1D coordinates (object)
[obj.Xo, obj.Yo] = meshgrid(obj.xo);% 2D coordinates (object)

% prevent negative indices bz centering spiral coordinates on object
R = R + obj.No/2 - obj.Np/2 + 50;
C = C + obj.No/2 - obj.Np/2 + 50;
obj.positions0 = [R',C'];

d = 1e-3;   % the smaller this parameter the larger the spatial frequencies in the simulated object
b = 47;     % topological charge (feel free to play with this number)
[theta, rho] = cart2pol(obj.Xo, obj.Yo);
t = ( 1 + sign( sin( b * theta + 2*pi * (rho/d).^2) ) )/2;
phaseFun = 1;
t = t .* circ(obj.Xo,obj.Yo, obj.Lo) .* (1 - circ(obj.Xo,obj.Yo,175*obj.dxo)) .* phaseFun + ...
    circ(obj.Xo,obj.Yo,150*obj.dxo) .* (1/4.*abs(1 + exp(1i*2*pi/(obj.wavelength*5e-2).*(obj.Xo.^2+obj.Yo.^2)) ).^2 > 1/2);
object = normconv2( t, gaussian2D(2,1) ); % smooth edges

% load object
figure(2); hsvxplot(object,'pixelSize',obj.dxo, 'colorbar', 'test')
axis image, title('object')

%% generate ptychogram (CP = diffraction data)

obj.ptychogram = zeros(obj.Nd, obj.Nd, obj.numFrames);
cmap = setColormap;
obj.propagator.type = 'Fresnel'; % specify propagator between sample and detector (Fraunhofer, Fresnel, ASP, scaledASP)
obj.propagator.quadraticPhase = exp(1i * 2*pi/obj.wavelength * (obj.Xp.^2 + obj.Yp.^2)/(2*obj.zo));
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
    if mod(loop,100) == 0
        figure(10)
        imagesc(obj.xd, obj.xd, log10(I + 1e-4)) 
        axis image; colormap(cmap); colorbar; 
        title('CP (diffraction) data')
        drawnow
    end
    
end

%%  calculate noise

% % simulate Poisson noise
% bitDepth = 16;
% maxNumCountsPerDiff = 2^bitDepth;
% 
% % normalize data (ptychogram)
% I = I/max(obj.ptychogram(:)) * maxNumCountsPerDiff;
% obj.ptychogram = obj.ptychogram/max(obj.ptychogram(:)) * maxNumCountsPerDiff;
% 
% % simulate Poisson noise
% % obj.ptychogram(:) = poisson_noise(obj.ptychogram(:));
% obj.ptychogram(:) = poissrnd(obj.ptychogram(:));
% 
% % compare noiseless data noisy 
% figure(12)
% imagesc(sqrt([I, obj.ptychogram(:,:,loop)])) 
% axis image off; colormap(cmap) 
% title(['left: noiseless, right: noisy (',num2str(bitDepth),' bit)'])

%% set data 

% note: the entrance pupil diameter determines the probe/pupil initial estimate diameter 
obj.entrancePupilDiameter = beamSize;

%% data inspection, check sampling requirements

obj.checkDataset

%% export CP data

obj.export.format = 'mat';
if exportBool
    obj.export.exportID = 'recentCP';
    obj.exportObj
end

%% convert to Fourier ptychography data
% I = reshape(permute(obj.ptychogram(:,:,1:end),[3,1,2]), [S, S, obj.Np^2]);
% I = flipud(rot90(I,1));
% 
% for k = round(obj.Np^2/2)+obj.Np/4:round(obj.Np^2/2) + obj.Np/2+10
% %    
%     figure(99)
%     imagesc(I(:,:,k))
%     axis image, colormap bone
%     title(num2str(k))
%     drawnow
% end

%%

% sin_phi_x = obj.Xd(:) ./ sqrt(obj.Xd(:).^2 + obj.Yd(:).^2 + obj.zo^2);
% sin_phi_y = obj.Yd(:) ./ sqrt(obj.Xd(:).^2 + obj.Yd(:).^2 + obj.zo^2);
% 
% obj.positions0 = -[sin_phi_y(:), sin_phi_x(:)]/obj.wavelength;


%%
obj.positions0 = zeros(obj.Np^2,2);
counter = 0;
for k = 1:obj.Np
    for l = 1:obj.Np
        counter = counter + 1;
        % compute angles 
        sin_phi_x = -obj.Xd(k,l) ./ sqrt(obj.Xd(k,l).^2 + obj.Yd(k,l).^2 + obj.zo^2);
        sin_phi_y = -obj.Yd(k,l) ./ sqrt(obj.Xd(k,l).^2 + obj.Yd(k,l).^2 + obj.zo^2);
        
        obj.positions0(counter,:) = [sin_phi_y, sin_phi_x];
    end
    
end
% convert angles to spatial frequencies
obj.positions0 = obj.positions0 / obj.wavelength;
%% convert to Fourier ptychography data and show result
obj.params.probe = imresize(fft2c(probe),[S,S]);
I = reshape(permute(obj.ptychogram(:,:,1:end),[3,1,2]), [S, S, obj.Np^2]);
% I = flipud(rot90(I,1));

for k = round(obj.Np^2/2)+obj.Np/4:round(obj.Np^2/2) + obj.Np/2+10    
    figure(99)
    imagesc(I(:,:,k))
    axis image, colormap bone
    title({'FP data', 'bright and dark field images'})
    drawnow
end

obj.operationMode = 'FPM'; % generate object (of class fracPty)
obj.ptychogram = I;
obj.numFrames = size(I,3);

%% this is the most important line of code in this tutorial:
% we can transform the CP data cube into a FP data cube with a single 
% line of code!

% obj.ptychogram = reshape(permute(obj.ptychogram(:,:,1:end),...
%                                     [3,1,2]), [S, S, obj.Np^2]);

%% exchange sampling

obj.Nd = S;
obj.dxd = dS * obj.dxp;         % detector pixel size
obj.Ld = obj.Nd * obj.dxd;      % detector size
obj.xd = (-obj.Nd/2:obj.Nd/2-1)*obj.dxd;% 1D coordinates (detector)
[obj.Xd, obj.Yd] = meshgrid(obj.xd);    % 2D coordinates (detector)
magnification = 1;                  % magnification
obj.zp = (magnification+1)*f;      	% pupil-to-detector = lens to image distance (optional, not mandatory)
obj.zo = obj.zp/magnification;      % object-to-lens distance (optional, not mandatory)

% object spectrum pixel size (reciprocal coordinates)
obj.dxo = 1/(obj.Ld/magnification); % object spectrum pixel size

% convert spatial frequencies to pixel units
obj.positions0 = round(obj.positions0 / obj.dxo);
% figure(99)
% scatter(obj.positions0(:,2), obj.positions0(:,1), 'ko', 'filled')
% xlabel('[px]'), ylabel('[px]'), axis square
% set(gca, 'FontSize', 20)
% axis square

%%

% entrance pupil (reciprocal coordinates)
obj.dxp = obj.dxo;              % pupil sampling step size
obj.Np = obj.Nd;                % number of samples in pupil
obj.Lp = obj.Np * obj.dxp;      % field of view in pupil
obj.xp = (-obj.Np/2:obj.Np/2-1)*obj.dxp;% 1D coordinates (pupil)
[obj.Xp, obj.Yp] = meshgrid(obj.xp);% 2D coordinates (pupil)

% object spectrum (reciprocal coordinates)
obj.No = 512;
obj.Lo = obj.No * obj.dxo;      % object spectrum field of view
obj.xo = (-obj.No/2:obj.No/2-0.5)*obj.dxo;% 1D coordinates (object spectrum)
[obj.Xo, obj.Yo] = meshgrid(obj.xo);    % 2D coordinates (object spectrum)

% center coordinates on object spectrum
obj.positions0 = obj.positions0 + obj.No/2 - obj.Np/2;
% add random offset
% obj.positions0 = obj.positions0 + round(2*(-0.5 + rand(size(obj.positions0))));

% figure(100)
% scatter(obj.positions0(:,2), obj.positions0(:,1), 'ko', 'filled')
% xlabel('[px]'), ylabel('[px]'), axis square
% set(gca, 'FontSize', 20)
% axis square


%% export FP data

obj.export.format = 'mat';
if exportBool
    obj.export.exportID = 'recentFP';
    obj.exportObj
end