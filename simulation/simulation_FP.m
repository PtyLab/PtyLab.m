% configuration

% clear
% close all

% note: set the paths in the config folder as needed (type 'open config' in command line)
s = pwd;
switch s(end-5:end)
    case 'ptyLab'
        addpath(genpath('config'))
    otherwise
        addpath(genpath('../config'))
end
[toolboxFolder, dataFolder] = configFile;   
addpath(genpath( toolboxFolder ));  % toolbox contains all required functions to run code
cd(toolboxFolder)

%% create object

obj = ptyLab('operationMode', 'FPM'); % generate object (of class ptyLab)
% note: type open 'ptyLab' to see the properties of this class
% - in the class definition you find required fills which are currently
% empty 
% - the purpose of this code is to fill in those properties and simulate 
% a CPM data set. To this end, all properties need to be defined. 

exportBool = true;              % if exportBool = true, data is saved in dataFolder
% get export path
obj.export.exportPath = dataFolder;    % this is where simulated data is saved

%% set physical properties

obj.wavelength = 632.8e-9;      % wavelength
obj.zo = 5e-2;                  % object pupil-distance distance

% detector coodinates (typically sampling requirements are set by the detector! So start from here)
obj.dxd = 2e-6;                 % detector pixel size
obj.Nd = 2^8;                   % number of samples in detector plane
obj.Ld = obj.Nd * obj.dxd;      % detector size
obj.xd = (-obj.Nd/2:obj.Nd/2-0.5)*obj.dxd;% 1D coordinates (detector)
[obj.Xd, obj.Yd] = meshgrid(obj.xd);    % 2D coordinates (detector)

% note: we assume a sampling grid equal in the source and object
% note: source coordinates are not part of the class
% note: we give the lens parameters first
f = 30e-3;              % focal length of singlet lens
M = 2;                  % magnification
obj.zp = (M+1)*f;      	% pupil-to-detector = lens to image distance
obj.zo = obj.zp/M;      % object-to-lens distance

% object spectrum pixel size (reciprocal coordinates)
obj.dxo = 1/(obj.Ld/M); % object spectrum pixel size

%% generate LED and Fourier-space positions

LEDspacing = 3e-3;  % LED spacing
z = 10e-2;          % LEd to object distance
n = 3;              % controls number of LEDs per dimension
disp(['number of LEDs: ', num2str((2*n+1)^2)])
[x_LED, y_LED] = meshgrid( (-n:n) * LEDspacing );
sin_phi_x = x_LED ./ sqrt(z^2 + x_LED.^2 + y_LED.^2); % direction sin -> x
sin_phi_y = y_LED ./ sqrt(z^2 + x_LED.^2 + y_LED.^2); % direction sin -> y

figure(1)
subplot(1,2,1)
scatter(sin_phi_x(:), sin_phi_y(:), 'ko', 'filled')
xlabel('[NA_x]'), ylabel('[NA_y]'), axis square
set(gca, 'FontSize', 20)

% use row-column order for Fourier space positions
obj.positions0 = [sin_phi_y(:), sin_phi_x(:)]/obj.wavelength;

% convert to pixel units and inflect
obj.positions0 = round(obj.positions0 / obj.dxo);
subplot(1,2,2)
scatter(obj.positions0(:,2), obj.positions0(:,1), 'ko', 'filled')
xlabel('[px]'), ylabel('[px]'), axis square
set(gca, 'FontSize', 20)
axis square

%% generate Fourier space coordinates

% entrance pupil (reciprocal coordinates)
obj.dxp = obj.dxo;              % pupil sampling step size
obj.Np = obj.Nd;                % number of samples in pupil
obj.Lp = obj.Np * obj.dxp;      % field of view in pupil
obj.xp = (-obj.Np/2:obj.Np/2-0.5)*obj.dxp;% 1D coordinates (pupil)
[obj.Xp, obj.Yp] = meshgrid(obj.xp);% 2D coordinates (pupil)

% object spectrum (reciprocal coordinates)
obj.No = max(obj.positions0(:)) - min(obj.positions0(:))+obj.Np+2;% number of samples in object spectrum plane
obj.No = 2^ceil(log2(obj.No));
obj.Lo = obj.No * obj.dxo;      % object spectrum field of view
obj.xo = (-obj.No/2:obj.No/2-0.5)*obj.dxo;% 1D coordinates (object spectrum)
[obj.Xo, obj.Yo] = meshgrid(obj.xo);    % 2D coordinates (object spectrum)

% center coordinates on object spectrum
obj.positions0 = obj.positions0 + obj.No/2 - obj.Np/2;
% add random offset
obj.positions0 = obj.positions0 + round(2*(-0.5 + rand(size(obj.positions0))));

figure(2)
scatter(obj.positions0(:,2), obj.positions0(:,1), 'ko', 'filled')
xlabel('[px]'), ylabel('[px]'), axis square
set(gca, 'FontSize', 20)
axis square

%% set source coordinates

dxs = 1/obj.Lo;        % source & object sampling
Ns = obj.No;              % number of source & object samples
Ls = Ns * dxs;          % source & object FOV
xs = (-Ns/2:Ns/2-0.5)*dxs;% 1D coordinates (source & object)
[Xs, Ys] = meshgrid(xs);    % 2D coordinates (source & object)


%% generate object

% option 1: spokes target
% d = 1e-4;   % the smaller this parameter the larger the spatial frequencies in the simulated object
% b = 13;     % topological charge (feel free to play with this number)
% [theta, rho] = cart2pol(Xs, Ys);
% t = ( 1 + sign( sin( b * theta + 2*pi * (rho/d).^2) ) )/2;
% phaseFun = exp(1i * theta);
% t = t .* circ(Xs,Ys, Ls/2) .* (1 - circ(Xs, Ys, Ls/16)) + ...
%     circ(Xs, Ys, Ls/32);
% object = convolve2( t, gaussian2D(5,1), 'same' ); % smooth edges

% option 2: cameraman
object = imresize( single( imread('cameraman.tif') ), obj.No * [1,1] );
phase = imresize( single( imread('westconcordorthophoto.png') ), obj.No * [1,1] );
phase = phase - mean(phase(:));
phase = phase / max(max(phase)) * pi;
object = object .* circ(obj.Xo, obj.Yo, obj.Lo/2) .* exp(1i * phase);

% load object
figure(3); hsvxplot(object,'pixelSize',xs, 'colorbar', 'test')
axis image, title('object')

%% generate pupil

obj.numFrames = length(obj.positions0);
pupilDiameter = 1/dxs/6/obj.dxo; % in pixel units (reused further below in those units)
aperture = circ(obj.Xo, obj.Yo, pupilDiameter * obj.dxo); 
H = circ(obj.Xo, obj.Yo, pupilDiameter * obj.dxo);
% note: the following line may be used to check the orientation of H in the
% reconstruction. If it's flipped, the data is incorrectly flipped.
% H = circ(obj.Xo, obj.Yo, pupilDiameter * obj.dxo) - circ(obj.Xo + pupilDiameter/3 * obj.dxo, obj.Yo, pupilDiameter/3 * obj.dxo);
% below some (strong) aberrations
jmax = 10;
[~, ~, Z] = ZernikeCalc( (1:jmax), H, H, 'Noll', [obj.No/2+1 obj.No/2+1 obj.No/2]);
H = H .* exp(1i * (0.5*Z(:,:,4) + 3*Z(:,:,5) + 0.5*Z(:,:,7))); % simulate aberrations
idx = obj.No/2 - obj.Np/2 +1:obj.No/2 + obj.Np/2;
cmap = setColormap;

figure(4); hsvxplot(H,'pixelSize',xs, 'colorbar', 'test')
axis image, title('object')

%% get data

for loop = 1:obj.numFrames
    r = sqrt(z^2 + (Xs - x_LED(loop)).^2 + (Ys - y_LED(loop)).^2);
%     r = sqrt(z^2 + (Xs + x_LED(loop)).^2 + (Ys + y_LED(loop)).^2);
    P = exp(1i * 2*pi/obj.wavelength * r);
    
    psi = object .* P;
    
    % get Fourier plane (high-NA) signal
    psi_tilde = fft2c(psi);
    
    % get low pass filter signal
    chi = fft2c( psi_tilde(idx, idx) .* H(idx,idx) );
    
    % get intensity
    obj.ptychogram(:,:,loop) = abs(chi).^2;
    
    figure(5)
    subplot(2,2,1)
    hsvxplot( psi(Ns/2-Ns/4:Ns/2+Ns/4, Ns/2-Ns/4:Ns/2+Ns/4) ,'pixelSize',dxs)
    title('exit wave (behind object)')
    
    subplot(2,2,2)
    imagesc( log(abs(psi_tilde)+1) ), axis image off, colormap(cmap)
    title('log-signal incident on pupil')
    
    subplot(2,2,3)
    hsvxplot(H)
    title('pupil')
    
    subplot(2,2,4)
    imagesc(obj.ptychogram(:,:,loop)), axis image off
    title('observed signal')
    drawnow
end

% estimate overlap

s = norm( obj.positions0(1,:) - obj.positions0(2,:), 2 );
overlap = 1 - s/pupilDiameter;
display(['linear overlap ', num2str(round(overlap * 1000)/10), ' %'])

%%  generate noise

% simulate Poisson noise
bitDepth = 8;
maxNumCountsPerDiff = 2^bitDepth;

idx = round(obj.numFrames/2);
I = round(obj.ptychogram(:,:,idx));
% normalize data (ptychogram)
I = I/max(obj.ptychogram(:)) * maxNumCountsPerDiff;
obj.ptychogram = obj.ptychogram/max(obj.ptychogram(:)) * maxNumCountsPerDiff;

% simulate Poisson noise
% obj.ptychogram(:) = poisson_noise(obj.ptychogram(:));
obj.ptychogram(:) = poissrnd(obj.ptychogram(:));

% compare noiseless data noisy 
figure(6)
imagesc(sqrt([I, obj.ptychogram(:,:,idx)])) 
axis image off; colormap(cmap) 
title(['left: noiseless,     right: noisy (',num2str(bitDepth),' bit)'])

%%  set entrance pupil

obj.entrancePupilDiameter = pupilDiameter * obj.dxo;

%% incflect data (this step would be done in preprocessing)

% obj.ptychogram = rot90(obj.ptychogram, 2);

%% set propagator
obj.propagator.type = 'Fraunhofer';

%% data inspection, check sampling requirements
obj.checkDataset

%% export data

obj.export.format = 'mat';
if exportBool
    obj.export.exportID = 'recent';
    obj.exportObj
end

%% TV test
% epsilon = 1e-2;s
% Gr = grad(H);
% d = sqrt(sum3(Gr.^2,3));
% G = -div( Gr ./ repmat( sqrt( epsilon^2 + d.^2 ) , [1 1 2]) );
% figure(99),hsvxplot(G)