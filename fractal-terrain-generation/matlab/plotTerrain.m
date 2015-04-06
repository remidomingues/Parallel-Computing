function plotTerrain(obj)
if ischar(obj)
    data = load(obj);
else
    data = obj;
end

figure
hold on
contourf(data,64,'EdgeColor','none','LineStyle','none')
surf(data,'EdgeColor','none','LineStyle','none','FaceLighting','phong');
colormap([seacolor;landcolor])
% axis tight
axis off
view([45 45]);
