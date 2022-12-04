clear
close all
tic
addpath(genpath('utils'))
%%


% rend = 200; ds = 30; Nr = round(rend / (ds)); [R, C] = GenerateConcentricCircles(Nr, ds, rend);
% ds = 15; [R, C] = GenerateArchimedeanSpiral(4, ds, 'c', 2);
[R, C] = GenerateFermatSpiral(150, 'c', 150);
% [R, C] = GenerateRasterGrid(10,1, 'type', 'randomOffset','amplitude', 0); R = 5 * R; C = 5 * C;
% [R, C] = GenerateRasterGrid(11, 10, 'type', 'randomOffset','amplitude', 0);
% R = [1 * ones(1,10), 2 * ones(1,10),3 * ones(1,10),4 * ones(1,10)];C = [(1:10), fliplr((1:10)), (1:10), fliplr((1:10))];

xy = [C', R'];
n = length(R); % number of points
distance = sqrt((R(1)-R(n))^2 + (C(1)-C(n))^2);
for k = 2:n
    distance = distance + sqrt((R(k)-R(k-1))^2 + (C(k)-C(k-1))^2);
end
display(['initial travel distance: ',num2str(distance)])

figure(99)
plot(C, R, 'ko-','MarkerFaceColor',[0,0,0]); hold on
plot(C(1), R(1), 'ro','MarkerFaceColor',[1 0 0]);
hold off
axis square

%%

pause(0.5)
userConfig = struct('xy',xy,'numIter',5e4, 'popSize', 40,'showResult',false); 
resultStruct = tsp_ga(userConfig);

%%

figure(100)
plot(C([1:n, 1]), R([1:n, 1]), 'o-','color',...
     0.5 * [1,1,1],'MarkerFaceColor',0*[1,1,1],'LineWidth',2); 
hold on
plot(C(1), R(1), 'ko','MarkerFaceColor',[1 0 0]);
axis square
% axis(5*[-105, 105, -105, 105])
title(['total distace: ', num2str( round(distance) )],'FontSize', 20)
hold off
set(gcf, 'Color', 'w');
% export_fig('originalGrid.pdf')

%%% 

optRoute = [resultStruct.optRoute, 1];
if length(optRoute) < 2000
    figure(101)
    plot(C(optRoute), R(optRoute), 'o--','color',...
        0.5 * [1,1,1],'MarkerFaceColor',0*[1,1,1],'LineWidth',2);
    hold on
    plot(C(1), R(1), 'ko','MarkerFaceColor',[1 0 0]);
    % title(['total distace: ', num2str( round(resultStruct.minDist) )],'FontSize', 20)
    axis square
    % axis(5*[-110, 110, -110, 110])
    hold off
    set(gcf, 'Color', 'w');
    % export_fig('optimizedGrid.pdf')
end

%% generate txt file

ColsForTXT = [R(optRoute); C(optRoute)];
averageDistance = sum(sqrt(diff(R(optRoute)).^2 + diff(C(optRoute)).^2)) / length(optRoute);
disp(['average step size: ', num2str(averageDistance)])
disp(['number of scan points: ', num2str(length(optRoute))])

% ColsForTXT = [R(1); C(1)];

fileID = fopen('positions.txt','w'); % open file for writing ('w')
fprintf(fileID, '%12s %4u\r\n', 'number of positions: ', size(ColsForTXT,1));
fprintf(fileID, '%12s %12s\r\n', 'y (row) [um] |','x (col) [um]');
fprintf(fileID, '%4.2f %4.2f\r\n', ColsForTXT);
fclose(fileID);

% x = 0:.1:1;
% A = [x; exp(x)];

% fileID = fopen('exp.txt','w');
% fprintf(fileID,'%6s %12s\n','x','exp(x)');
% fprintf(fileID,'%6.2f %12.8f\n',A);
% fclose(fileID);

toc

