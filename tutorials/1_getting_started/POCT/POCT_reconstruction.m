% configuration
set(0,'DefaultFigureWindowStyle','normal')
clear
close all
% when array size is prime, use fftw planner
% fftw('planner','patient');

%% add toolbox and data folder

directory_info = dir();
toolboxFolder = erase( directory_info(1).folder, 'reconstruction');
dataFolder  = 'D:\Du\ptyLabExport';

addpath(genpath( toolboxFolder ));
addpath(genpath( dataFolder ));
cd(toolboxFolder)

%% load data set

load('TwoLayer_bin4')

%% initial params

% initial estimate of beam size
obj.entrancePupilDiameter = obj.Np/3 * obj.dxp; 
obj = obj.initialParams(...
    'initialProbe', 'circ',...
    'initialObject', 'ones',...
    'saveMemory', false,...
    'deleteFrames', [ ],...
    'fontSize', 17);

%% manual params
obj.params.objectPlot = 'abs';              % 'complex', 'piComlex', 'abs', 'angle', or 'piAngle' 
obj.params.intensityConstraint = 'Anscombe';% 'Anscombe', 'ProxAnscombe', 'Poisson', 'ProxPoisson', 'Gaussian', 'BoseEinstein'
obj.params.numIterations = 1000;

% main parameters
obj.params.betaObject = 0.25;               % gradient step size object
obj.params.betaProbe = 0.25;                % gradient step size probe
obj.params.npsm = 1;                        % number of orthogonal modes
obj.params.FourierMaskSwitch = false;       % apply mask to corrupted pixels
obj.params.gpuSwitch = true;                % gpuSwitch

% object regularization 
obj.params.objectSmoothenessSwitch = true;  % if true, impose smootheness
obj.params.objectSmoothenessWidth = 2;      % # pixels over which object is assumed fairly smooth
obj.params.objectSmoothnessAleph = 1e-2;    % relaxation constant that determines strength of regularization

% probe regularization 
obj.params.probePowerCorrectionSwitch = true;   % probe normalization to measured PSD
obj.params.comStabilizationSwitch = true;       % center of mass stabilization for probe

% other parameters
obj.propagator.type = 'Fraunhofer';     % specify propagator between sample and detector (Fraunhofer, Fresnel, ASP, scaledASP)
obj.params.verboseLevel = 'high';       % control how much output is plotted
obj.params.fftshiftSwitch = true;       % if true, the FFT is not centered
obj.monitor.objectZoom = 1.1;             % allows for object zoom in
obj.monitor.objectMax = 1;              % readjust relative intensity

% mPIE reconstruction
obj.monitor.figureUpdateFrequency = 1; % frequency of reconstruction monitor 
obj.params.engine = 'mPIE';             % 'ePIE', 'mPIE', 'zPIE, 'e3PIE', 'm3PIE', 'pcPIE', 
obj = reconstruct(obj);                 % run reconstruction

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
end
