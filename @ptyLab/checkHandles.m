function checkHandles(obj)
% this function initializes handles 

figure(2)
subplot(1,2,1)
set(gcf, 'Color', 'w');
cmap = colormap(obj.monitor.DIFFcmap);
maxI = log10( max(obj.ptychogram(:)) + 1);
if obj.params.fftshiftSwitch
    obj.monitor.fig2_1 = imagesc( log10( fftshift(obj.ptychogram(:,:,1)) + 1 ), [0, maxI] );
else
    obj.monitor.fig2_1 = imagesc( log10( gather(obj.ptychogram(:,:,1)) + 1 ), [0, maxI] );
end
colormap(cmap)
h = colorbar;
ylabel(h, 'log_{10}(counts)')
obj.monitor.ax2_1 = gca;
xlabel('[px]'), ylabel('[px]')
axis image
set(gca,'FontSize',obj.monitor.fontSize)


subplot(1,2,2)
if obj.params.CPSCswitch
    cScale = size(obj.ptychogram, 1) / size(obj.probe, 1);
    if obj.params.fftshiftSwitch
        obj.monitor.fig2_2 = imagesc( obj.xd, obj.xd, log10( fftshift( obj.ptychogram(:,:,1)) + 1 ), [0, maxI + log10(cScale)] );
    else 
        obj.monitor.fig2_2 = imagesc( obj.xd, obj.xd, log10( obj.ptychogram(:,:,1) + 1 ), [0, maxI + log10(cScale)] );
    end
else
    if obj.params.fftshiftSwitch
        obj.monitor.fig2_2 = imagesc(obj.xd, obj.xd, fftshift(log10( obj.ptychogram(:,:,1) )), [0, maxI] );
    else
        obj.monitor.fig2_2 = imagesc(obj.xd, obj.xd, obj.ptychogram(:,:,1), [0, maxI] );
    end
end
h = colorbar;
ylabel(h, 'log_{10}(counts)')
title({['I_m @ ', num2str( round(obj.zo *1000) ), ' mm'],''})
axis image
set(gca,'FontSize',obj.monitor.fontSize)