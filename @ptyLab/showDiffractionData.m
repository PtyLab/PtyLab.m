function showDiffractionData(obj)
% showDiffractionData:      show estimated and measured diffraction patterns

if obj.params.fftshiftSwitch
    obj.monitor.fig2_1.CData = gather( log10( fftshift(obj.params.Iestimated) + 1 ) );
else
    obj.monitor.fig2_1.CData = log10( gather(obj.params.Iestimated) + 1 );
end

obj.monitor.ax2_1.Title.String = {['I_e @ #: ', num2str(obj.params.currentPosition)],''};

if obj.params.CPSCswitch
    if obj.params.fftshiftSwitch
        obj.monitor.fig2_2.CData = gather( log10( fftshift(obj.params.ImeasuredDownsampled) ) );
    else 
        obj.monitor.fig2_2.CData = gather( log10( obj.params.ImeasuredDownsampled ) );
    end
else
    if obj.params.fftshiftSwitch
        obj.monitor.fig2_2.CData = gather( log10( fftshift(obj.params.Imeasured) ) );
    else
        obj.monitor.fig2_2.CData = gather( log10(obj.params.Imeasured) );
    end
end

end