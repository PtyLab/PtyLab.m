% note: this script recontructs the freshly exported FP data set from the
% simulation in the same folder

% configuration
set(0,'DefaultFigureWindowStyle','normal')
clear
close all
restoredefaultpath
set(0, 'DefaultAxesFontSize', 18);
set(0, 'DefaultFigureColor', 'w');
set(0,'defaultAxesFontName', 'serif')

%% add toolbox and data folder

directory_info = dir();
toolboxFolder = erase( directory_info(1).folder, 'tutorials\5_conversion_CP_FP');
dataFolder  = 'C:\Users\Lars Loetgering\Documents\ptyLabExport';

addpath(genpath( toolboxFolder ));
addpath(genpath( dataFolder ));
cd(toolboxFolder)

%% load data set

load('recentFP')

% subtract a bit of background
% obj.ptychogram = nonnegative( obj.ptychogram - 100 );

%% initialize
% set pupil radius (= detection NA)
% obj.params.NA = obj.Lp/4 * obj.wavelength;
% obj.entrancePupilDiameter = obj.params.NA * obj.dxo;

obj = obj.initialParams(...
    'initialProbe', 'circ',...
    'initialObject', 'ones',...
    'saveMemory', false,...
    'deleteFrames', [ ],...
    'fontSize', 17);

%% manual initialization

obj.object = ifft2c(imresize(sqrt( mean(obj.ptychogram, 3) ), ...
                    obj.No * [1,1]) );
obj.probe = circ(obj.Xp, obj.Yp, obj.Lp/2);

%% manual params

obj.params.objectPlot = 'complex';          % 'complex', 'abs', 'angle', 'piAngle'
obj.params.intensityConstraint = 'Anscombe';% 'Anscombe', 'ProxAnscombe', 'Poisson', 'ProxPoisson', 'Gaussian', 'BoseEinstein'

% engine
obj.params.engine = 'pSD';  % parallel steepest descent

% main parameters
obj.monitor.figureUpdateFrequency = 1;  % frequency of reconstruction monitor 
obj.params.numIterations = 10000;       % total number of iterations
obj.params.betaObject = 1;              % gradient step size object
obj.params.betaProbe = 1;               % gradient step size probe
obj.params.npsm = 1;                    % number of probe/pupil state mixtures
obj.params.FourierMaskSwitch = false;   % apply mask to corrupted pixels
obj.params.gpuSwitch = true;            % gpuSwitch

% object regularization 
obj.params.objectSmoothenessSwitch = true;      % if true, impose smootheness
obj.params.objectSmoothenessWidth = 3;          % # pixels over which object is assumed fairly smooth
obj.params.objectSmoothnessAleph = 1e-2;        % relaxation constant that determines strength of regularization
obj.params.absObjectSwitch = false;             % force the object to be abs-only
obj.params.absObjectBeta = 1e-2;                % relaxation parameter for abs-only constraint
obj.params.objectContrastSwitch = false;        % pushes object to zero outside ROI
obj.params.objectTVregSwitch = false;           % TV regularization on object

% probe regularization 
obj.params.probeSmoothenessSwitch = true;       % enforce probe smootheness (prevent grid artefact) 
obj.params.probeSmoothnessAleph = 1e-2;         % relaxation parameter for probe smootheness
obj.params.probeSmoothenessWidth = 3;           % loose object support diameter
obj.params.absorbingProbeBoundary = false;      % controls if probe has period boundary conditions (zero)
obj.params.probePowerCorrectionSwitch = false;  % probe normalization to measured PSD
obj.params.modulusEnforcedProbeSwitch = false;  % enforce empty beam
obj.params.comStabilizationSwitch = false;      % center of mass stabilization for probe
obj.params.probeTVregSwitch = false;

% other parameters
obj.params.positionOrder = 'sequential';% 'sequential' or 'random', 'FP'
obj.propagator.type = 'Fraunhofer'; % specify propagator between sample and detector (Fraunhofer, Fresnel, ASP, scaledASP)
obj.params.backgroundModeSwitch = false;    % background estimate
obj.params.makeGIF = false;         % export GIF animation of reconstruction
obj.params.orthogonalizationFrequency = 100; % # iterations until orthogonalization
obj.params.noprpsm = 4;             % number of modes for for (i)OPRP
obj.params.verboseLevel = 'high';   % control how much output is plotted
obj.params.fftshiftSwitch = true;   % if true, no  time is lost on fftshifting
obj.params.probeROI = [1 obj.Np 1 obj.Np];  % probe region of interest
obj.monitor.objectZoom = 1;                 % control object zoom in
obj.monitor.objectMax = 1;                  % control maximum brightness in object reconstruction
obj.monitor.showObjectSpectrum = false;     % show object spectrum
obj = reconstruct(obj);
