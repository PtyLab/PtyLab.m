function showPtychogram(obj, varargin)
% show ptychogram on a log10 - scale

p = inputParser;
p.parse( varargin{:} )
cmap = setColormap;

if isfield(obj.params, 'fftshiftSwitch')
    if obj.params.fftshiftSwitch
        sliceViewer( fftshift(log10(obj.ptychogram + 1)), 'Colormap', cmap );
    else
        sliceViewer( log10(obj.ptychogram + 1), 'Colormap', cmap );
    end
else
    sliceViewer( log10(obj.ptychogram + 1), 'Colormap', cmap );
end

title('ptychogram')

screen_size = get(0,'screensize');
h = gcf;
h.Position(1:2) = round( min(screen_size(3:4)))/4*[1,1];
h.Position(3:4) = round( min(screen_size(3:4)))/2*[1,1];
end
