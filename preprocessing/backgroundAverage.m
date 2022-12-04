% generate background average
clear
restoredefaultpath
% add toolbox folders
switch getenv('COMPUTERNAME')
    case 'PHASESPACE'
        toolboxFolder = 'F:\Dropbox\Codes\fracPty'; addpath(genpath(toolboxFolder)); 
        obj = fracPty; obj.export.exportPath = 'F:\fracPtyExport';
        
    case and('', isunix)
        toolboxFolder = '/home/lars/ownCloud/IXO/fracPty'; addpath(genpath(toolboxFolder)); 
        obj = fracPty; obj.export.exportPath = '/home/lars/Documents/fracPtyExport';
    case 'EUVLAB5'
        toolboxFolder = '\\storage01\data\ARCNL\groups\eikema-witte-group\Phasespace\ptychography\fracPty'; addpath(genpath(toolboxFolder)); 
        obj = fracPty; obj.export.exportPath = '\\storage01\data\ARCNL\groups\eikema-witte-group\Phasespace\ptychography\fracPty\datasets\ptychoData';
    otherwise
        error('computer and paths not properly specified')
end

cd(toolboxFolder);
% change this depending on computer used for preprocessing
filePath = '\\sun.amolf.nl\eikema-witte\group-folder\Phasespace\ptychography\rawData\USAF_dorian\20200506_162352_VocationalGuidance00016\AVT camera (GT3400)';

N = 1456; % number of detector pixels

%%

% number of Frames
temp = dir(filePath);
numFrames = length(temp) - 2;

backgroundAv = 0;

disp('==== read data ====')
for k = 0:numFrames-1
    if k<10
        fileNameTemp = ['0000',num2str(k),'.tif'];
    elseif and( 9 < k, k < 100)
        fileNameTemp = ['000',num2str(k),'.tif'];
    elseif and( 99 < k, k < 1000)
        fileNameTemp = ['00',num2str(k),'.tif'];
    elseif and( 999 < k, k < 10000)
        fileNameTemp = ['0',num2str(k),'.tif'];
    end
    disp(['...reading frame ',fileNameTemp])
    
    backgroundAv = backgroundAv + posit(single(imread([filePath,'/',fileNameTemp])));
end

backgroundAv = backgroundAv/numFrames;

%%

figure(1)
subplot(1,2,1)

imagesc(posit(single(imread([filePath,'/',fileNameTemp]))))
axis image
colorbar
colormap gray
title('left: single frame')

subplot(1,2,2)
imagesc(backgroundAv)
axis image
colorbar
colormap gray
title('right: average background')

%%

imwrite( uint16(backgroundAv), [filePath,'/background.tif'])