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

% first argument: number of points per dimension
% second argument: step size in um
n = 30;
ds = 80;
[R, C] = GenerateRasterGrid(n, ds, 'amplitude', 10);

% undo flipping

for k = 2:2:n
    
    counter = (k-1) * n + 1;
    R(counter:counter + n-1) = fliplr(R(counter:counter + n-1));
    C(counter:counter + n-1) = fliplr(C(counter:counter + n-1));
end


% R = [0, R, 0];
% C = [0, C, 0];
idx = abs(C) < 750;
R = [0, R(idx), 0];
C = [0, C(idx), 0];
% R = [R,R(1)];
% C = [C,C(1)];

figure(99)
plot(C, R, 'ko-','MarkerFaceColor',[0,0,0]); hold on
plot(C(1), R(1), 'ro','MarkerFaceColor',[1 0 0]);
hold off
axis equal
set(gcf,'color','w')

optRoute = 1:length(R); % number of points
averageDistance = sum(sqrt(diff(R(optRoute(2:end-1))).^2 + diff(C(optRoute(2:end-1))).^2)) / length(optRoute(2:end-1));
disp(['averageDistance: ',num2str(averageDistance)])

%% generate txt file
ColsForTXT = [R(optRoute); C(optRoute)];
disp(['number of scan points: ', num2str(length(optRoute))])

distance = 0;
for k = 2:n
    distance = distance + sqrt((R(optRoute(k))-R(optRoute(k-1)))^2 + (C(optRoute(k))-C(optRoute(k-1)))^2);
end
disp(['final travel distance: ',num2str(distance)])

fileID = fopen('positions.txt','w'); % open file for writing ('w')
fprintf(fileID, '%12s %4u\r\n', 'number of positions: ', size(ColsForTXT,1));
fprintf(fileID, '%12s %12s\r\n', 'y (row) [um] |','x (col) [um]');
fprintf(fileID, '%4.2f %4.2f\r\n', ColsForTXT);
fclose(fileID);

toc

