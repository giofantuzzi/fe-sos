function plotHermite(x,u,N,varargin)

% Plot hermite interpolation with N points for each element
nel = length(x)-1;
uplot = zeros(nel*N,1);
xplot = zeros(nel*N,1);

% Interpolate on all elements
shift = 0;
pos = 1:N;
idx = 1:4;
for i = 1:nel
    % Mesh on element
    a = x(i);
    b = x(i+1);
    xloc = linspace(a,b,N).';
    xplot(pos) = xloc;
    % Basis functions
    bfun =  [-((b - xloc).^2.*(b - 3*a + 2.*xloc))./(a - b)^3, ...
             -((a - xloc).*(b - xloc).^2)./(a - b)^2, ...
             ((a - xloc).^2.*(a - 3*b + 2.*xloc))/(a - b)^3, ...
             -((a - xloc).^2.*(b - xloc))./(a - b)^2];
    % Interpolate
    uplot(pos) = bfun*u(idx);
    % Shift
    pos = pos + N;
    idx = idx + 2;
end

% Unique and plot
[xplot,i] = unique(xplot);
h = ishold;
pp = plot(xplot,uplot(i),varargin{:}); hold on
% plot(x,u(1:2:end),'.','Color',pp.Color);
if ~h; hold off; end
