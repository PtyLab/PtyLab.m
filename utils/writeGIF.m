function writeGIF(filename, write_bool, counter, delayTime, varargin)

p = inputParser;
p.addParameter('figNum', 1)
p.parse( varargin{:} );

if write_bool
    frame = getframe(p.Results.figNum);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im, 256);
    if counter == 1
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime', delayTime);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime', delayTime);
    end
end
end