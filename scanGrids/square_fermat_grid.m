% configuration
set(0, 'DefaultAxesFontSize', 18);
set(0, 'DefaultFigureColor', 'w');
clc
clear
close all
restoredefaultpath
tic
% add relevant folders
switch getenv('COMPUTERNAME')
    case 'DESKTOP-MCS5K1H'
        toolboxFolder = 'C:\Users\Lars Loetgering\Dropbox\Codes\ptyLab';
end
        
addpath(genpath( toolboxFolder ));
cd(toolboxFolder);

% set true to export results
export_bool = true;

%% generate (non-optimal) grid

% first argument: number of points, 
% second argument: scaling parameter
[R, C] = generate_fermat_spiral(800, 'scaling', 20);
[R, C] = rectify_fermat(R,C);

travel_distance = sum(sqrt(diff(R).^2 + diff(C).^2));
disp(['initial travel distance: ',num2str(round(travel_distance))])

figure(1)
plot(C, R, 'ko-','MarkerFaceColor',[0,0,0]); hold on
plot(C(1), R(1), 'ro','MarkerFaceColor',[1 0 0]);
hold off
axis square tight off

%% solve traveling salesman problem via 2opt

[Ropt, Copt] = twoOpt(R,C);

travel_distance_opt = sum(sqrt(diff(Ropt).^2 + diff(Copt).^2));
disp(['initial travel distance: ',num2str(round(travel_distance_opt))])

%% show original grid

cmap = turbo(length(R));
sz = 50;

figure(1)
plot(C(1:end-1), R(1:end-1), 'k-','MarkerFaceColor',[0,0,0]); hold on
scatter(C(1:end-1), R(1:end-1), sz, cmap(1:end-1,:), 'filled')
hold off
axis square tight off

% show optimization result

cmap = turbo(length(Ropt));

figure(2)
plot(Copt(1:end-1), Ropt(1:end-1), 'k-','MarkerFaceColor',[0,0,0]); hold on
scatter(Copt(1:end-1), Ropt(1:end-1), sz, cmap(1:end-1,:), 'filled')
hold off
axis square tight off

%% generate txt file
if export_bool
    ColsForTXT = [Ropt; Copt];
    averageDistance = sum(sqrt(diff(Ropt).^2 + diff(Copt).^2)) / length(Ropt);
    disp(['final travel distance: ',num2str(round(averageDistance * length(Ropt)))])
    disp(['average step size [um]: ', num2str(averageDistance)])
    disp(['number of scan points: ', num2str(length(Ropt))])
    % open file for writing ('w')
    fileID = fopen([toolboxFolder,'\scanGrids\generated_grids\positions.txt'],'w'); 
    fprintf(fileID, '%12s %4u\r\n', 'number of positions: ', size(ColsForTXT,1));
    fprintf(fileID, '%12s %12s\r\n', 'y (row) [um] |','x (col) [um]');
    fprintf(fileID, '%4.2f %4.2f\r\n', ColsForTXT);
    fclose(fileID);
end

