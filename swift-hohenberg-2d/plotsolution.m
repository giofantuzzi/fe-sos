function plotsolution(Xp,Yp,Up,TRI,x,y,plotMesh)

% Plot solution at interpolated points with mesh overlaid
if nargin < 7; plotMesh = 1; end
L = max(max(Yp));
C = 2.25;%ceil(max(max(Up)));
plt = surf(Xp,Yp,Up);
plt.LineStyle = 'none';
view([0 90]);
if plotMesh
    hold on
    t = triplot(TRI,x,y);
    t.LineWidth = 0.2;
    t.Color = [1 1 1]*0.9;
    t.ZData = 3*ones(size(t.XData));
    hold off
end
axis equal
xlim([-2*L,2*L])
ylim([-L,L])
colormap(jet)
caxis([-1,2])
% colorbar
% drawnow
end