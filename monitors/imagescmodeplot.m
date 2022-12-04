function imagescmodeplot(P, varargin)

M = size(P,1);  % number of rows
N = size(P,2);  % number of columns in stack
s = floor(sqrt(size(P,3))); % number of slices shown (must be quadratic)

% parse optional inputs
p = inputParser;
normalizeFlag = true;
p.addParameter('normalize', normalizeFlag)
p.addParameter('pixelSize', 1)
p.addParameter('contrastModality', 'abs')
p.parse(varargin{:})

Q = zeros(s * size(P,1), s * size(P,2), 'like', P);
counter = 0;
for k = 1:s
    for l = 1:s
        counter = counter + 1;
        if p.Results.normalize
            Ptemp = P(:,:,counter);
            Q((k-1)*M+1:(k-1)*M+M, (l-1)*N+1:(l-1)*N+N) = Ptemp / (max(abs(Ptemp(:)))+eps);
        else
            Q((k-1)*M+1:(k-1)*M+M, (l-1)*N+1:(l-1)*N+N) = P(:,:,counter)  ;
        end
    end
end

m = size(Q,1);  % number of rows
n = size(Q,2);  % number of columns in stack
x = ((1:n)-round(n/2)) * p.Results.pixelSize;
y = ((1:m)-round(m/2)) * p.Results.pixelSize;

switch p.Results.contrastModality
    case 'angle'
        temp = angle(Q);
        imagesc(x,y,temp - min(temp(:)))
        colormap bone
    case 'abs'
        imagesc(x,y,abs(Q),[0,1])
        colormap bone
        h = colorbar('Location','West');
        h.FontSize = 20;
    case 'piAngle'
        temp = angle(Q);
        n = 1;
        imagesc(x,y,temp - mean(temp(:)),[-pi/n,pi/n])
        colormap(gray)
        h = colorbar('Location','West');
        h.Ticks = [-pi/n,0,pi/n];
        h.TickLabels = {num2str(round(-pi/n*10)/10); '\phi=0'; ['+',num2str(round(pi/n*10)/10)]};
        h.Color = [1 1 1];
        
%         colormap(bluewhitered(200))
%         colormap(othercolor('Accent8'))
%         colormap(othercolor('Paired10'))
%         colormap(othercolor('BuGy_8'))
%         imagesc(x,y,temp - mean(temp(:)),[-pi/2,pi/2])
        
    otherwise
        error('contrast modality for imagescmodeplot incorrectly specified')
end


axis image
set(gcf, 'Color', 'w');
return
