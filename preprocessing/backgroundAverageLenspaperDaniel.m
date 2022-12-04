% generate background average
clear
restoredefaultpath
% add toolbox folders
switch getenv('COMPUTERNAME')
    case and('', isunix)
        % TODO: change this according to your the local path containing the
        % background measurements
        toolboxFolder = '/home/user/Dropbox/Codes/ptyLab'; addpath(genpath(toolboxFolder)); 
     otherwise
        error('computer and paths not properly specified')
end

cd(toolboxFolder);
% change this depending on computer used for preprocessing
filePath = '/home/user/Dropbox/CurrentProjects/ptychoWorkshop/data/27-05-20/bk';

N = 1040; % number of detector pixels
% N = 1456; % number of detector pixels

%%

% number of Frames
temp = dir(filePath);
numFrames = length(temp) - 3;

backgroundAv = 0;

disp('==== read data ====')
for k = 0:numFrames-1
    fileNameTemp = ['bk',num2str(k+1),'.png'];
    disp(['...reading frame ',fileNameTemp])
    temp = posit(single(rgb2gray(imread([filePath,'/',fileNameTemp]))));
    backgroundAv = backgroundAv + temp;
end

backgroundAv = backgroundAv/numFrames;

%%

figure(1)
subplot(1,2,1)

imagesc(temp)
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