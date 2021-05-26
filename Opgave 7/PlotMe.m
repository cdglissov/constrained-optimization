function PlotMe(xlimits,ylimits)

x = linspace(xlimits(1), xlimits(2),100);
y = linspace(ylimits(1), ylimits(2),100);
[X,Y] = meshgrid(x,y);

F  = (X.^2+Y-11).^2 + (X+Y.^2-7).^2;
v = [0:1:10 10:10:100 100:10:200];

figure(1);
[c,h]=contour(X,Y,F,v,'linewidth',1);
colorbar
axis image

% Constraints with y isolated
yc1 = (x+2).^2;
yc2 = 4/10*x;

hold on
xlim([-5 ,5])
ylim([-5, 5])
fill(x,yc1,[0.01 0.01 0.01],'facealpha',0.2)
plot(x,yc2,'linewidth',0.5, 'color' , [0.001 0.001 0.001]);
plot(x,yc1,'linewidth',0.5, 'color' , [0.001 0.001 0.001]);
fill([x x(end) x(1)],[yc2 -5 -5],[0.01 0.01 0.01],'facealpha',0.2)
