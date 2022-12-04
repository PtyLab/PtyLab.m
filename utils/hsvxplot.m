function varargout = hsvxplot(u, varargin)
% hsvplot(u) generates hue-brightness plot of two dimensional input u
% last change: 21st August 2018

p = inputParser;
p.addParameter('intensityScale', [ ])
p.addParameter('colorbar', false);
p.addParameter('scalebar', 0);
p.addParameter('pixelsize', 0)
p.parse( varargin{:} );

% normalize birghtness (value) to range [0, 1]
u = gather(u);
r = abs(u);
r = r / ( max(r(:)) + eps );
if ~isempty(p.Results.intensityScale)
    r = posit(r - p.Results.intensityScale(1));
    r = r/(p.Results.intensityScale(2) - p.Results.intensityScale(1));
end

% normalize angle
phi = angle( u );
phi = ( phi + pi )/( 2 * pi );

% normalization of phase saturation
B = zeros(size(u,1), size(u,2), 3);         %Declare RGB array
B(:,:,1) = phi;
B(:,:,2) = 1;
B(:,:,3) = r;
A = hsv2rgb(B);

switch p.Results.colorbar
    
    case 'both'
        % draws both hue and brightness bars (left and bottom)
        k = ceil(size(u,1)/15);
        
        % phase bar
        m = size(u,1);
        D = zeros(m, 1, 3);
        D(:,:,1) = (0:1/m:1-1/m)';
        D(:,:,2) = 1;
        D(:,:,3) = 1;
        C = hsv2rgb(D);
        C(:,1:k,:) = repmat( C(:,1,:),1, k );
        A=[C, A];
        
        % intensity bar
        n = size(A, 2);
        E = zeros(1, n, 3);
        E(:,:,1) = 1;
        E(:,:,2) = 0;
        E(:,:,3) = (0:1/n:1-1/n);
        F = hsv2rgb(E);
        F = repmat( F, k, 1 );
        A = [A; F];
        
    case 'North'
        % draws hue bar north 
        m = size(u,1);
        k = ceil(m/15);
        
        % phase bar
        D = zeros(1, m, 3);
        D(:,:,1) = (0:1/m:1-1/m);
        D(:,:,2) = 1;
        D(:,:,3) = 1;
        C = hsv2rgb(D);
        C = repmat( C, k, 1 );
        A(1:k, :, :) = C;
        
    case 'East'
        % draws hue bar east
        m = size(u,1);
        k = ceil(m/15);
        
        % phase bar
        D = zeros(m, 1, 3);
        D(:,:,1) = (0:1/m:1-1/m)';
        D(:,:,2) = 1;
        D(:,:,3) = 1;
        C = hsv2rgb(D);
        C = repmat( C, 1, k );
        A(:, m-k+1:m, :) = C;
        
    case 'South'
        % draws hue bar south
        m = size(u,1);
        k = ceil(m/15);
        
        % phase bar
        D = zeros(1, m, 3);
        D(:,:,1) = (0:1/m:1-1/m);
        D(:,:,2) = 1;
        D(:,:,3) = 1;
        C = hsv2rgb(D);
        C = repmat( C, k, 1 );
        A(m-k+1:m, :, :) = C;
        
    case 'West'
        % draws hue bar west
%         m = size(u,1);
%         k = ceil(m/15);
%         
%         % phase bar
%         D = zeros(m, 1, 3);
%         D(:,:,1) = (0:1/m:1-1/m)';
%         D(:,:,2) = 1;
%         D(:,:,3) = 1;
%         C = hsv2rgb(D);
%         C = repmat( C, 1, k );
%         A(:, 1:k, :) = C;

        m = size(u,1);
        k = ceil(m/15);
        
        % phase bar
        D = zeros(m, 1, 3);
        D(:,:,1) = (0:1/m:1-1/m)';
        D(:,:,2) = 1;
        D(:,:,3) = 1;
        C = hsv2rgb(D);
        
    case 'test'
        m = size(u,1);
        k = ceil(m/15);
        
        % phase bar
        D = zeros(m, 1, 3);
        D(:,:,1) = (0:1/m:1-1/m)';
        D(:,:,2) = 1;
        D(:,:,3) = 1;
        C = hsv2rgb(D);
        
end

if p.Results.scalebar > 0
    
    m = round(size(A, 1)/20);
    n = round(size(A, 2)/20);
    scaleBar = ones( m + 1, p.Results.scalebar + 1, 3 ,'like', A);
    scaleBar(:,:,1) = 0.1;
    A(m:2*m,(3*n):(3*n+p.Results.scalebar),:) = hsv2rgb(scaleBar);
    
end

if p.Results.pixelsize > 0
    N = size(A,2); x = (1:N) * p.Results.pixelsize;
    M = size(A,1); y = (1:M) * p.Results.pixelsize;
    imagesc(x,y,A); axis square
else
    imagesc(A); axis image off;
end
set(gcf, 'Color', 'w');
%
if strcmp(p.Results.colorbar,'test')
    colormap(squeeze(C))
%     c = colorbar('Location','East');
%     set(c,'YTick',[])
    c = colorbar('Location','East','Ticks',[0.05,0.95],'TickLabels', {'\color{white} 0', '\color{white} 2\pi'});
    set(gca,'FontSize',17)
%     c.TickLabels = {'0', '2\pi'} ;
elseif strcmp(p.Results.colorbar,'West')
    colormap(squeeze(C))
    c = colorbar('Location','West','Ticks',[0.05,0.95],'TickLabels', {'\color{white} 0', '\color{white} 2\pi'});
    c.TickLabels = [];
    % 	set(gca,'FontSize',20)
end
%

switch nargout
    case 1
        varargout{1}=uint8(A*255);
    otherwise
end

return