function plotTerrain(varargin)
plotType = 0;
data = [];

for i=1:nargin
    if ischar(varargin{i})
        if (strcmp(varargin{i},'3d') || strcmp(varargin{i},'surf'))
            plotType = 1;
        elseif (strcmp(varargin{i},'flat') || strcmp(varargin{i},'image'))
            plotType = 2;
        else
            if isempty(data)
                data = load(varargin{i});
            else
                data = data + load(varargin{i});
            end
        end
    else
        if isempty(data)
            data = varargin{i};
        else
            data = data + varargin{i};
        end
    end
end

if (plotType == 1 || plotType == 0)
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
    axis off;
    axis equal;
end
end
