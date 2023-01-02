% The goal of this tutorial is to show zPIE in action.
% 
% Maximize the reconstruction window to monitor how the object gradually
% appears better as the object-detector distance converges.
% 
% You can play with "obj.params.zPIEgradientStepSize" to get a feeling for
% the interplay between stability versus convergence rate of zPIE. 

% configuration
set(0,'DefaultFigureWindowStyle','normal')
clear
close all
% when array size is prime, use fftw planner
% fftw('planner','patient');

%% add toolbox and data folder

directory_info = dir();
toolboxFolder = erase( directory_info(1).folder, 'reconstruction');
dataFolder  = 'C:\Users\Lars Loetgering\Documents\ptyLabExport';

addpath(genpath( toolboxFolder ));
addpath(genpath( dataFolder ));
cd(toolboxFolder)

%% load data set

load('recent')

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
obj.params.numIterations = 200;

% main parameters
obj.params.betaObject = 0.75;               % gradient step size object
obj.params.betaProbe = 0.75;                % gradient step size probe
obj.params.npsm = 2;                        % number of orthogonal modes
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
obj.monitor.objectZoom = 2;             % allows for object zoom in
obj.monitor.objectMax = 1;              % readjust relative intensity

% mPIE reconstruction
obj.monitor.figureUpdateFrequency = 10; % frequency of reconstruction monitor 
obj.params.engine = 'zPIE';             % 'ePIE', 'mPIE', 'zPIE, 'e3PIE', 'm3PIE', 'pcPIE', 
obj.params.zPIEgradientStepSize = 200;  % step size for zPIE
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
