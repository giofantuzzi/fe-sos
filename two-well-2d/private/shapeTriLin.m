function [phi,gradxiphi]=shapeTriLin(xieta)
%Special case of linear triangle has constant gradient
xi=xieta(:,1);
eta=xieta(:,2);
phi(:,1)=xi;
phi(:,2)=eta;
phi(:,3)=1-xi-eta;
gradxiphi(:,1)=[1;0];
gradxiphi(:,2)=[0;1];
gradxiphi(:,3)=[-1;-1];