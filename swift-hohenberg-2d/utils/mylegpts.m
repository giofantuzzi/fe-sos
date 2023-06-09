function [x,w] = mylegpts(n)
% Calculates Legendre roots and weights for Gauss-Legendre quadrature using
% the traditional Golub-Welsch eigenvalue method (less effective than most
% other alternatives, but of historical value).
beta=.5./sqrt(1-(2*(1:n-1)).^(-2)); %3-term recurrence coeffs
T=diag(beta,1)+diag(beta,-1);       %Jacobi matrix
[V,D]=eig(T);                       %Eigenvalue decomposition
x=diag(D);                          %Legendre points
[x,i]=sort(x);                      %Sort
w=2*V(1,i).^2;                      %Quadrature weights
%Enforce symmetry:
ii=1:floor(n/2);
x=x(ii);
w=w(ii);
if (mod(n,2))%Odd number.
    x=[x;0;-x(end:-1:1)];
    w=[w,2-sum(2*w),w(end:-1:1)];
else %Even number.
    x=[x;-x(end:-1:1)];
    w=[w,w(end:-1:1)];
end
end