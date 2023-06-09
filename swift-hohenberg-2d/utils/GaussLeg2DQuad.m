function [xi2D,w2D]=GaussLeg2DQuad(n)
[xi,w1D]=mylegpts(n);
[Xi,Eta]=meshgrid(xi,xi);
W2D=w1D'*w1D;
xi2D=[Xi(:),Eta(:)];
w2D=W2D(:)';
