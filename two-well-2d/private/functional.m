function [E,S,N]=EnergyLinElem(u,nodes,elements,constants)
% Compute "energy"
% \int_Omega 0.5*c*|grad(u)|^2 + (u-a)^2*(u-b)^2 dOmega
IEN=elements';
n_el=size(IEN,2); %Number of elements (also n_el=size(elements,1))
a=constants(1); b=constants(2); c=constants(3);
E=0;
S=0;
N=0;
%Loop over elements
for i=1:n_el
    xnod=nodes(IEN(:,i),:);
    alphaloc=u(IEN(:,i));
    %First the 0.5*c*|grad(u)|^2 using nP=2 and element stiffness matrix
    ke=localElemStiffLin(xnod);
    S=S+0.5*c*alphaloc'*ke*alphaloc;
    %Then the higher order nonlinear term
    nonlin=@(z) ((z-a).^2).*((z-b).^2);
    nP=4; %Order of nonlinearity
    [xi2D,w2D]=PrecomputedGaussLeg2DTri(nP);
    [phi,gradxiphi]=shapeTriLin(xi2D);
    [~,detJ,~]=shapeTriElemLin(xnod,phi,gradxiphi);
    uloc=phi*alphaloc;
    N=N+detJ*w2D*nonlin(uloc);
end
E = N+S;