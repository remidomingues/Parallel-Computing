function plotTerrain(obj)
if ischar(obj)
    data = load(obj);
else
    data = obj;
end

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
