% configuration
set(0,'DefaultFigureWindowStyle','normal')

clear
close all
restoredefaultpath
tic

% add relevant folders
switch getenv('COMPUTERNAME')
    case 'PHASESPACE'
        toolboxFolder = 'F:\Dropbox\Codes\fracPtyPoly';
        dataFolder = 'F:\fracPtyExport';
    case and('', isunix)
        toolboxFolder = '/home/lars/Dropbox/Codes/fracPtyPoly';
        dataFolder = '/home/lars/Documents/fracPtyExport';
    case 'GOYA'
        toolboxFolder = 'C:\Users\PR\Documents\MATLAB\fracPtyPoly';addpath(genpath(toolboxFolder));
        obj = fracPty; obj.exportPath = 'C:\Users\PR\Documents\MATLAB\fracPtyExport';
end
        
addpath(genpath( toolboxFolder ));
addpath(genpath( dataFolder ));
cd(toolboxFolder);

%% parameters

numPoints = 300;        % number of points
radius = 22;            % radius of final scan grid [in micrometers]
p = 1.25;               % "clumping parameter" (see notes* below)
beamSize = 10e-6;       % expected beam size, required to calculate overlap
numIterations = 1e5;    % number of iterations in optimization

% * p = 1 is standard Fermat 
% * p > 1 yields more points towards the center of grid   

%% generate (non-optimal) grid

[R, C] = GenerateFermatSpiral(numPoints);
n = length(R); % number of points
distance = sqrt((R(1)-R(n))^2 + (C(1)-C(n))^2);
for k = 2:n
    distance = distance + sqrt((R(k)-R(k-1))^2 + (C(k)-C(k-1))^2);
end
disp(['initial travel distance: ',num2str(distance)])

%%

maxD = max( sqrt(R.^2 + C.^2) );
R = R/maxD; 
C = C/maxD; 
r = sqrt(R.^2 + C.^2);
theta = atan2(R, C);
R = radius * r.^p .* cos(theta);
C = radius * r.^p .* sin(theta);
 
figure(1)
h = gcf;
h.Color = [1,1,1];
h.Units = 'normalized';
h.OuterPosition = [1 0 1 1];
subplot(2,2,1)
plot(C, R, 'ko-','MarkerFaceColor', [0,0,0]); hold on
plot(C(1), R(1), 'ro','MarkerFaceColor', [1 0 0]);
hold off
axis square
title('initial scan grid')
h = gca;
h.XLim = radius * [-1,1];
h.YLim = radius * [-1,1];
xy = [C', R'];

%% traveling salesman problem, genetic algorithm (tsp_ga)

userConfig = struct('xy',xy,'numIter', numIterations, 'popSize', 50); 
disp('optimizing scan grid ...')
resultStruct = tsp_ga(userConfig);

figure(1),subplot(2,2,2)
h = gca;
h.XLim = radius * [-1,1];
h.YLim = radius * [-1,1];

%% generate txt file

optRoute = [resultStruct.optRoute, 1];
ColsForTXT = [R(optRoute); C(optRoute)];
distances = sqrt(diff(R(optRoute)).^2 + diff(C(optRoute)).^2);
averageDistance = sum(distances) / length(optRoute);
disp(['average step size [um]: ', num2str(averageDistance)])
disp(['number of scan points: ', num2str(length(optRoute))])

%% overlap

overlap = 1-distances*1e-6/beamSize;
figure(1),subplot(2,2,3)
histogram(100*overlap, 30)
axis square, grid on
xlabel('overlap [%]')
ylabel('absoslute frequency')
h = gca;
h.XLim = [50 100];

figure(1),subplot(2,2,4)
for k = 1:numPoints
    cmap = (overlap(k) > 0.8) * [0 1 0] + ...
        and(overlap(k) > 0.7, overlap(k) < 0.8) * [0 0 1] + ...
        (overlap(k) < 0.7) * [1 0 0];
    plot(R(optRoute(k)), C(optRoute(k)),'o','MarkerFaceColor',cmap,'MarkerEdgeColor',cmap)
    hold on
    axis square
end
hold off
title({'green: o > 80%', 'blue: 70% < o < 80%', 'red: o < 70%'})
h = gca;
h.XLim = radius * [-1,1];
h.YLim = radius * [-1,1];

%%

fileID = fopen('positions.txt','w'); % open file for writing ('w')
fprintf(fileID, '%12s %4u\r\n', 'number of positions: ', size(ColsForTXT,1));
fprintf(fileID, '%12s %12s\r\n', 'y (row) [um] |','x (col) [um]');
fprintf(fileID, '%4.2f %4.2f\r\n', ColsForTXT);
fclose(fileID);
toc

