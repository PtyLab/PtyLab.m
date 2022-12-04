function r = setColormap

cmap = [[0,0,0]; jet(20); [1,1,1]];

n = 100;
x = linspace(0,1,size(cmap, 1));
xq = linspace(0,1,n);

cmap2(:,1) = interp1(x,cmap(:,1),xq,'linear');
cmap2(:,2) = interp1(x,cmap(:,2),xq,'linear');
cmap2(:,3) = interp1(x,cmap(:,3),xq,'linear');

r = cmap2;

end