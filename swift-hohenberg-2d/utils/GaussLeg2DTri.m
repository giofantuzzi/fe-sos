function [xi2D,w2D]=GaussLeg2DTri(nP)
[xi2DQ,w2DQ]=GaussLeg2DQuad(nP);
xi2D=[0.5*(xi2DQ(:,1)+1).*(1-0.5*(xi2DQ(:,2)+1)),0.5*(xi2DQ(:,2)+1)];
w2D=0.25*w2DQ.*(1-0.5*(xi2DQ(:,2)'+1));