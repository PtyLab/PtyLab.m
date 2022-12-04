function varargout = hsvmodeplot(P, varargin)

 % parse optional inputs
 p = inputParser;
 p.addParameter('normalize', true);
 p.addParameter('pixelSize',0);
 p.parse(varargin{:});
 
if size(P,3) > 3
    M = size(P,1);  % number of rows
    N = size(P,2);  % number of columns in stack
    s = floor(sqrt(size(P,3))); % number of slices shown (must be quadratic)
    
    Q = zeros(s * size(P,1), s * size(P,2),'like',P);
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
    
    if p.Results.pixelSize > 0
        A = hsvplot(Q);
        M = size(A, 1); N = size(A, 2);
        y = (1:M) * p.Results.pixelSize;
        x = (1:N) * p.Results.pixelSize;  
        imagesc(x,y,A); axis image
    else
        hsvplot(Q)
    end
    
else
    hsvxplot(P(:,:,1),'pixelsize',p.Results.pixelSize)
end

if nargout==1
    varargout{1} = Q;
end
return
