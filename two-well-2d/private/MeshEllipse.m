function [nodes,elements,dirichlet]=ElliptMesh(h)
% Mesh ellipse using triangles of size h
a=1; b=0.5; %major and minor axes of ellipse
hr=a/ceil(a/h); %mesh size (<=h)
[X,Y]=meshgrid(-a:hr:a,-b:hr:b);
X=X(:); Y=Y(:); XY=[X,Y];
%Remove nodes outside the domain: 0.5*hr away from boundary
tol=4*eps;
XYIC=XY(((XY(:,1)/a).^2+(XY(:,2)/b).^2<(1-0.5*hr+tol)^2),:);
%Boundary discretization
thetah=2*pi/ceil(2*pi*a/(sqrt(2)*h)); %theta spacing
theta=0:thetah:(round(2*pi/thetah)-1)*thetah;
%Final set of nodes: interior + boundary
nodes=[XYIC;a*cos(theta(:)),b*sin(theta(:))];
elements=delaunay(nodes(:,1),nodes(:,2)); %convex hull meshing
nodalbdryind=size(XYIC,1)+1:size(nodes,1); %only last nodes are bdry
dirichlet=[nodalbdryind(:),zeros(length(nodalbdryind),1)];
end