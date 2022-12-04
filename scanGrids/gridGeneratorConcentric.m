% configuration
set(0,'DefaultFigureWindowStyle','normal')

clear
close all
restoredefaultpath
tic
% add relevant folders
switch getenv('COMPUTERNAME')
    case 'PHASESPACE'
        toolboxFolder = 'F:\Dropbox\Codes\fracPty';
        dataFolder = 'F:\fracPtyExport';
    case and('', isunix)
        toolboxFolder = '/home/lars/Dropbox/Codes/fracPty';
        dataFolder = '/home/lars/Documents/fracPtyExport';
    case 'EUVSTUDENT1'
        obj.exportPath = 'D:\fracPtyExport';
end
        
addpath(genpath( toolboxFolder ));
addpath(genpath( dataFolder ));
cd(toolboxFolder);

%% generate (non-optimal) grid

% first argument: number of points, second argument: scaling of Fermat grid
Nr = 40;
s = 70;
rend = 3500;

% Nr = 20;
% s = 70;
% rend = 1000;
[R, C] = GenerateConcentricGrid(Nr, s, rend);

xy = [C', R'];
n = length(R); % number of points
distance = sqrt((R(1)-R(n))^2 + (C(1)-C(n))^2);
for k = 2:n
    distance = distance + sqrt((R(k)-R(k-1))^2 + (C(k)-C(k-1))^2);
end
disp(['initial travel distance: ',num2str(distance)])

figure(99)
plot(C, R, 'ko-','MarkerFaceColor',[0,0,0]); hold on
plot(C(1), R(1), 'ro','MarkerFaceColor',[1 0 0]);
hold off
axis square

%%

rng('default'); % make random offset reproducible
randomAmplitude = 5;
C = C + round( randomAmplitude * ( -1 + 2 * rand(size(C))));
R = R + round( randomAmplitude * ( -1 + 2 * rand(size(R))));



C(1) = 0; C(end) = 0;
R(1) = 0; R(end) = 0;

%% show optimization result

figure(100)
plot(C([1:n, 1]), R([1:n, 1]), 'o-','color',...
     0.5 * [1,1,1],'MarkerFaceColor',0*[1,1,1],'LineWidth',2); 
hold on
plot(C(1), R(1), 'ko','MarkerFaceColor',[1 0 0]);
axis((rend+randomAmplitude)*[-1 1 -1 1]), axis square
title(['total distace: ', num2str( round(distance) )],'FontSize', 20)
hold off
set(gcf, 'Color', 'w');

%% generate txt file

ColsForTXT = [R; C];
averageDistance = sum(sqrt(diff(R).^2 + diff(C).^2)) / length(C);
disp(['average step size [um]: ', num2str(averageDistance)])
disp(['number of scan points: ', num2str(length(C))])

fileID = fopen('positions.txt','w'); % open file for writing ('w')
fprintf(fileID, '%12s %4u\r\n', 'number of positions: ', size(ColsForTXT,1));
fprintf(fileID, '%12s %12s\r\n', 'y (row) [um] |','x (col) [um]');
fprintf(fileID, '%4.2f %4.2f\r\n', ColsForTXT);
fclose(fileID);

toc
