function [xy,detJ,gradxphi]=shapeTriElemLin(xnod,phi,gradxiphi)
%Special case of linear triangle has constant gradient elementwise
xy=zeros(size(phi,1),2);
J=[0,0;0,0];
for c=1:3
    xy=xy+kron(phi(:,c),xnod(c,:));
    J=J+xnod(c,:)'*gradxiphi(:,c)';
end
detJ=det(J);
gradxphi=inv(J).'*gradxiphi;