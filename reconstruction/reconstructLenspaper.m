% configuration
set(0,'DefaultFigureWindowStyle','normal')

clear
close all
restoredefaultpath

% add relevant folders
s = pwd;

switch s(end-5:end)
    case 'ptyLab'
        addpath(genpath('config'))
    otherwise
        addpath(genpath('../config'))
end

[toolboxFolder, dataFolder] = configFile;
addpath(genpath( toolboxFolder ));
addpath(genpath( dataFolder ));
cd(toolboxFolder)

%% load data set

load('recent')

%%

obj.entrancePupilDiameter = obj.Np/3 * obj.dxp; % initial estimate of beam
obj = obj.initialParams(...
    'initialProbe', 'circ',...
    'initialObject', 'ones',...
    'saveMemory', false,...
    'deleteFrames', [ ],...
    'fontSize', 17);

obj.probe = obj.probe .* exp(1i * 2*pi/obj.wavelength * (obj.Xp.^2 + obj.Yp.^2)/(2*6e-3));

figure(1)
hsvplot(obj.probe)

%% manual params

obj.params.objectPlot = 'complex';          % 'complex', 'piComlex', 'abs', 'angle', or 'piAngle'
obj.params.intensityConstraint = 'standard';% 'standard', 'exponential', fluctuation

% engine
obj.params.engine = 'mPIE';                 % 'ePIE', 'mPIE', 'zPIE, 'e3PIE', 'm3PIE', 'pcPIE', 'kPIE', 'OPRP', 'SD', 'msPIE', 'sDR'

% main parameters
obj.monitor.figureUpdateFrequency = 5;   % frequency of reconstruction monitor
obj.params.numIterations = 1000;        % total number of iterations
obj.params.betaObject = 0.25;           % gradient step size object
obj.params.betaProbe = 0.25;            % gradient step size probe
obj.params.npsm = 4;                    % number of orthogonal modes
obj.params.FourierMaskSwitch = false;   % apply mask to corrupted pixels
obj.params.gpuSwitch = false;            % gpuSwitch

% object regularization
obj.params.objectSmoothenessSwitch = false;     % if true, impose smootheness
obj.params.objectSmoothenessWidth = 2;          % # pixels over which object is assumed fairly smooth
obj.params.objectSmoothnessAleph = 1e-2;        % relaxation constant that determines strength of regularization
obj.params.absObjectSwitch = false;             % force the object to be abs-only
obj.params.absObjectBeta = 1e-2;                % relaxation parameter for abs-only constraint
obj.params.objectContrastSwitch = true;        % pushes object to zero outside ROI

% probe regularization
obj.params.probeSmoothenessSwitch = false;      % enforce probe smootheness
obj.params.probeSmoothnessAleph = 5e-2;         % relaxation parameter for probe smootheness
obj.params.probeSmoothenessWidth = 3;           % loose object support diameter
obj.params.absorbingProbeBoundary = false;      % controls if probe has period boundary conditions (zero)
obj.params.probePowerCorrectionSwitch = true;   % probe normalization to measured PSD
obj.params.modulusEnforcedProbeSwitch = false;  % enforce empty beam
obj.params.comStabilizationSwitch = true;       % center of mass stabilization for probe

% other parameters
obj.params.positionOrder = 'random';    % 'sequential' or 'random'
obj.propagator.type = 'Fraunhofer';     % specify propagator between sample and detector (Fraunhofer, Fresnel, ASP, scaledASP)
obj.params.backgroundModeSwitch = false;    % background estimate
obj.params.makeGIF = false;         % export GIF animation of reconstruction
obj.params.orthogonalizationFrequency = 10; % # iterations until orthogonalization
obj.params.noprpsm = 4;             % number of modes for for (i)OPRP
obj.params.verboseLevel = 'high';   % control how much output is plotted
obj.params.fftshiftSwitch = true;
obj.params.probeROI = [1 obj.Np 1 obj.Np];
obj.monitor.objectZoom = 1;
obj.monitor.objectMax = 0.5;
obj.monitor.showObjectSpectrum = false;
obj = reconstruct(obj);

% Pcalibrated = gather(obj.probe(:,:,1)); save('Pcalibrated','Pcalibrated')
% Ocalibrated = gather(obj.object); save('Ocalibrated','Ocalibrated')

%% export data

% switch to CPU
obj.params.gpuSwitch = false;
obj.checkGPU

% export data
exportBool = false;
if exportBool
    % set path for data export
    obj.getExportPath;
    % set file name for data small export
    obj.exportID = [obj.fileName,'_recSmall'];
    % export essential data
    obj.exportReconstruction
    % set file name for data large export
    obj.exportID = [obj.fileName,'_recLarge'];
    % export data
    obj.exportObj
    % exit to release licenses (relevant for institutions with limited licenses)
%     exit
end
