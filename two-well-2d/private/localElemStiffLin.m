function ke=localElemStiffLin(xnod)
%With linear elements, gradients are constant elementwise, so it is only
%necessary to multiply by area of integration (detJ*0.5, where 0.5 is area
%of master element).
[xi2D,w2D]=PrecomputedGaussLeg2DTri(1);%Argument is irrelevant (cnst grad)
[phi,gradxiphi]=shapeTriLin(xi2D);
[~,detJ,gradxphi]=shapeTriElemLin(xnod,phi,gradxiphi);
ke=gradxphi'*gradxphi*detJ*sum(w2D); 