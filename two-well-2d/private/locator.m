%It computes the ID,IEN and LM arrays as described by Hughes..
function [ID,IEN,LM]=locator(nodes,elements,dirichlet)
Nn=size(nodes,1);
ND=size(dirichlet,1);
ff=ones(1,Nn);
ff(dirichlet(:,1))=0;
ff=logical(ff);
ID=zeros(1,Nn);
ID(ff)=1:(Nn-ND);
IEN=transpose(elements);
LM=ID(IEN);