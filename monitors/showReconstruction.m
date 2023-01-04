function showReconstruction(obj)
% showReconstruction:       main plot for reconstruction

%% object plot

figure(1)
subplot(1,3,1)
    
switch obj.operationMode
    case 'CPM'
        object = obj.object(obj.params.objectROI(1):obj.params.objectROI(2), obj.params.objectROI(3):obj.params.objectROI(4));
        dxo = obj.dxo;
    case 'FPM'
        object = ifft2c( center(obj.object) );
        dxo = 1/obj.Lo;
end

switch obj.params.objectPlot
    
    case 'complex'
        complex_plot(object,'colorbar', 'inside', 'pixelsize', dxo,'intensityScale',[0, obj.monitor.objectMax])
    case 'abs'
        imagescmodeplot( object, 'pixelSize', dxo )
    case 'angle'
        imagescmodeplot( object, 'pixelSize', dxo, 'contrastModality', 'angle' ), colorbar('Location','West')
    case 'piAngle'
        imagescmodeplot( object, 'pixelSize', dxo, 'contrastModality', 'piAngle' )
        
end

axis image
zoom(obj.monitor.objectZoom)
title({'object',' '}, 'FontWeight','normal')
set(gca,'FontSize',obj.monitor.fontSize)

%% probe plot %%

subplot(1,3,2)
hsvmodeplot(obj.probe(obj.params.probeROI(1):obj.params.probeROI(2),obj.params.probeROI(3):obj.params.probeROI(4),:),...
    'normalize', true,'pixelSize',obj.dxp)

% remove quadratic phase
% probeTemp = bsxfun(@times, obj.probe, conj(obj.ODpropagator.quadraticPhase));
% hsvmodeplot(probeTemp(obj.params.probeROI(1):obj.params.probeROI(2),obj.params.probeROI(3):obj.params.probeROI(4),:),...
%     'normalize', true,'pixelSize',obj.dxp)

switch obj.operationMode
    case 'CPM'
        ptitle = 'probe';
    otherwise
        ptitle = 'pupil';
end
            

if isfield(obj.params,'purity')
    % initialize error
    title({ptitle,['purity: ',num2str(round(obj.params.purity*100)),' %']},'FontWeight','normal')
else
    title({ptitle,' '},'FontWeight','normal')
end
set(gca,'FontSize',obj.monitor.fontSize)

%% error plot %%

if length(obj.params.error) >= 2^4
    idx = [2.^(1:floor(log2(length(obj.params.error)))), length(obj.params.error)];
else
    idx = 1:length(obj.params.error);
end
subplot(1,3,3)
loglog(idx, obj.params.error(idx),'-o')
axis square; grid on
title({'error',' '},'FontWeight','normal')
set(gca,'FontSize',obj.monitor.fontSize)
drawnow

%% other plots

% additional FP plot
if obj.monitor.showObjectSpectrum
    figure(3)
    switch obj.operationMode
        case 'CPM'
            temp = fft2c( obj.object(obj.params.objectROI(1):obj.params.objectROI(2), ...
                                     obj.params.objectROI(3):obj.params.objectROI(4), 1 ) );
            imagesc(log(abs( temp ) + 1)); axis image, colormap gray
        case 'FPM'
            imagesc(log(abs(obj.object(obj.params.objectROI(1):obj.params.objectROI(2), ...
                obj.params.objectROI(3):obj.params.objectROI(4),: ))+1)); axis image, colormap gray
    end
    drawnow
end

if obj.params.backgroundModeSwitch && strcmp(obj.params.verboseLevel, 'high')
    figure(4)
    set(gcf,'color','w')
    if obj.params.fftshiftSwitch
        imagesc(log10( real(fftshift(obj.params.background)) + 1 ))
        colormap(obj.monitor.DIFFcmap)
    else
        imagesc(log10( real(gather(obj.params.background)) + 1 ))
        
    end
    axis image, colormap(obj.monitor.DIFFcmap), colorbar
    title('background')
    drawnow
end

disp('===========================')
disp(['iteration:',num2str(length(obj.params.error))])
disp(['runtime: ',num2str(toc),' seconds'])
disp(['error: ',num2str(obj.params.error(end),'%10.3e\n')])
if isfield(obj.params,'purity')
    s = [];
    %     for k = 1:min(obj.params.npsm, 9)
    for k = 1:min(length(obj.params.normalizedEigenvaluesProbe), 25)
        s = [s, num2str(round(1000*obj.params.normalizedEigenvaluesProbe(k))/10),' '];
    end
    disp(['coherence structure: ', s, ' [%] '])
end

if strcmp(obj.params.verboseLevel, 'high')
    
%     if ~obj.params.saveMemory
%         obj.showErrorAtPosition
%     end
    obj.showDiffractionData
    obj.getOverlap(1, 2)
end

if obj.params.makeGIF
    if ~isfield(obj.params,'GIFcounter')
        obj.params.GIFcounter = 1;
    else
        obj.params.GIFcounter = obj.params.GIFcounter + 1;
    end
    
    if ~isfield(obj.params,'GIFfilename')
        obj.params.GIFfilename = 'recentGIF';
    end
    
    writeGIF(obj.params.GIFfilename, obj.params.makeGIF, obj.params.GIFcounter, 0.1)
end

return