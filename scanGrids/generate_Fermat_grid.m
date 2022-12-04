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
[R, C] = GenerateFermatSpiral2(150, 'c', 30);

xy = [C', R'];
n = length(R); % number of points
distance = sqrt((R(1)-R(n))^2 + (C(1)-C(n))^2);
for k = 2:n
    distance = distance + sqrt((R(k)-R(k-1))^2 + (C(k)-C(k-1))^2);
end
disp(['initial travel distance: ',num2str(round(distance))])

figure(99)
plot(C, R, 'ko-','MarkerFaceColor',[0,0,0]); hold on
plot(C(1), R(1), 'ro','MarkerFaceColor',[1 0 0]);
hold off
axis square

%% traveling salesman problem, genetic algorithm (tsp_ga)

pause(0.5)
userConfig = struct('xy',xy,'numIter',1e5, 'popSize', 40, 'showProg', true); 
resultStruct = tsp_ga(userConfig);

%% show optimization result

figure(100)
plot(C([1:n, 1]), R([1:n, 1]), 'o-','color',...
     0.5 * [1,1,1],'MarkerFaceColor',0*[1,1,1],'LineWidth',2); 
hold on
plot(C(1), R(1), 'ko','MarkerFaceColor',[1 0 0]);
axis square
title(['total distace: ', num2str( round(distance) )],'FontSize', 20)
hold off
set(gcf, 'Color', 'w');

%%% 

optRoute = [resultStruct.optRoute, 1];
figure(101)
set(gcf, 'Color', 'w');
if length(optRoute) < 2000
    
    plot(C(optRoute), R(optRoute), 'o-','color',...
        0.5 * [1,1,1],'MarkerFaceColor',0*[1,1,1],'LineWidth',2);
    hold on
    plot(C(1), R(1), 'ko','MarkerFaceColor',[1 0 0]);
    axis square
    hold off
    
end
grid on
xlabel('\mum')
ylabel('\mum')
%% generate txt file

ColsForTXT = [R(optRoute); C(optRoute)];
averageDistance = sum(sqrt(diff(R(optRoute)).^2 + diff(C(optRoute)).^2)) / length(optRoute);
disp(['final travel distance: ',num2str(round(averageDistance * length(optRoute)))])
disp(['average step size [um]: ', num2str(averageDistance)])
disp(['number of scan points: ', num2str(length(optRoute))])


fileID = fopen('positions.txt','w'); % open file for writing ('w')
fprintf(fileID, '%12s %4u\r\n', 'number of positions: ', size(ColsForTXT,1));
fprintf(fileID, '%12s %12s\r\n', 'y (row) [um] |','x (col) [um]');
fprintf(fileID, '%4.2f %4.2f\r\n', ColsForTXT);
fclose(fileID);

toc

