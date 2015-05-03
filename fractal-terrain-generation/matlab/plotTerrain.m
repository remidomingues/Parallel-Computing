function plotTerrain(varargin)
if ischar(varargin{1})
    data = load(varargin{1});
else
    data = varargin{1};
end

plotType = 0;
if nargin > 1,
    if (strcmp(varargin{2},'3d') || strcmp(varargin{2},'surf'))
        plotType = 1;
    elseif (strcmp(varargin{2},'flat') || strcmp(varargin{2},'image'))
        plotType = 2;
    else fprintf('Unknown option %s\n', varargin{2});
    end
end

if plotType <= 1
    n = 256;
    hmax = max(max(data));
    hmin = min(min(data));
    
    nland = round(hmax*n/(hmax-hmin))+1;
    nsea = n - nland;
    figure
    hold on
    contourf(data,32,'EdgeColor','none','LineStyle','none')
    surf(data,'EdgeColor','none','LineStyle','none','FaceLighting','phong');
    colormap([seacolor(nsea);landcolor(nland)])
    axis off
    view([45 45]);
    
elseif (plotType == 2 || plotType == 0)
    imagesc(data);
    colormap(copper);
    colorbar;
end
end
