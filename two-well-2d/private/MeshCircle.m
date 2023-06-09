function [nodes,elements,dirichlet]=CircleMesh(h)
% Mesh unit circle using triangles of size h
R=1; %radius
hr=R/ceil(R/h); %mesh size (<=h)
[X,Y]=meshgrid(-R:hr:R);
X=X(:); Y=Y(:); XY=[X,Y];
%Remove nodes outside the domain: 0.5*hr away from boundary
tol=4*eps;
XYIC=XY((XY(:,1).^2+XY(:,2).^2<(R-0.5*hr+tol)^2),:);
%Boundary discretization
thetah=2*pi/ceil(2*pi*R/(sqrt(2)*h)); %theta spacing
theta=0:thetah:(round(2*pi/thetah)-1)*thetah;
%Final set of nodes: interior + boundary
nodes=[XYIC;R*cos(theta(:)),R*sin(theta(:))];
elements=delaunay(nodes(:,1),nodes(:,2)); %convex hull meshing
nodalbdryind=size(XYIC,1)+1:size(nodes,1); %only last nodes are bdry
dirichlet=[nodalbdryind(:),zeros(length(nodalbdryind),1)];
end