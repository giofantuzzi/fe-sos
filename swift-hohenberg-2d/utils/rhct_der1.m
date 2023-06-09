function K = rhct_der1(X,C,d)

% Calculate "stiffness" matrix \int Du' * C * Du dx on a reduced HCT
% element defined with the 3 points in the cols of X. C is given. NOTE:
% Setting C = [1 0; 0 1] corresponds to the norm of the gradient since we use 
% the convetion that Du = (u_x, u_y)
%
% TAKEN FROM: https://www.degruyter.com/view/journals/cmam/12/4/article-p486.xml

% ---------------------------------------------------------------------------- %
% Gauss quadrature
% ---------------------------------------------------------------------------- %
% cub2D = gaussData(d);
% gaussPoints = 0.5.*( cub2D(:,1:2) + 1 );
% gaussWeights = cub2D(:,3)./4;
[gaussPoints,gaussWeights] = PrecomputedGaussLeg2DTri(d);
persistent Dxi dold
if isempty(Dxi) || d~=dold
    dold = d;
    for i = 1:length(gaussWeights)
       Dxi{i} = evalbasis(gaussPoints(i,:));
    end
end
% ---------------------------------------------------------------------------- %
% PRELIMINARIES
% ---------------------------------------------------------------------------- %
% Triangle edges
E = X(:,[3 1 2]) - X(:,[2 3 1]);

% Normal to edges
R = [0 -1; 1 0];
N = R*E;

% inter-element edges
F = ( E(:,[2 3 1]) - E(:,[3 1 2]) )./3;

% Volume
J = [E(:,2), E(:,3)];
mu = abs(det(J))/3;

% vectors b{k,j}
b = cell(3,3);
for k = 1:3
    [kp,km] = rotateindex(k);
    normEk2 = dot(E(:,k),E(:,k));
    boo = 3*mu/normEk2;
    b{k,kp} = [6*dot(E(:,k), F(:,km))/normEk2; 2*F(:,km)+boo*N(:,k)];
    b{k,km} = [-6*dot(E(:,k), F(:,kp))/normEk2; 2*F(:,kp)+boo*N(:,k)];
end

% Matrices S and T
S = 6.*[3, 3, 3; F].';
M = cell(3,1);
for j = 1:3
    [jp,jm] = rotateindex(j);
    T = zeros(3,3);
    T(jm,:) = b{jp,j}.';
    T(jp,:) = b{jm,j}.';
    T(j,:) = (b{jp,j} + b{jm,j} + [6; -2*F(:,j)]).';
    M{j} = S\T;
end

% ---------------------------------------------------------------------------- %
% ASSEMBLE
% ---------------------------------------------------------------------------- %
K = zeros(9,9);
D{1} = zeros(2,3);
D{2} = zeros(2,3);
D{3} = zeros(2,3);
for k = 1:3
   [kp,km] = rotateindex(k);
   J = F(:,[kp,km]);
   H = [1, 0 0; [0;0], J.'];
   G = [J(2,2), -J(2,1); -J(1,2), J(1,1)]./mu;
   for i = 1:length(gaussWeights)
       %Dxi = evalbasis(gaussPoints(i,:));
       tmp = Dxi{i}{1}*H;
       D{k}  = tmp*M{k};
       D{kp} = tmp*M{kp} + Dxi{i}{2}*H + Dxi{i}{4}*b{k,kp}.';
       D{km} = tmp*M{km} + Dxi{i}{3}*H + Dxi{i}{4}*b{k,km}.';
       D2 = G*[D{1}, D{2}, D{3}];
       K = K + (gaussWeights(i)*mu) .* ( D2.' * C * D2);
   end
end
K = 0.5*(K+K.');
K(abs(K)<1e-12) = 0;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [kp,km] = rotateindex(k)
kp = k+1;
km = k-1;
if k==1; km=3; end
if k==3; kp=1; end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Dxi = evalbasis(x)
y = x(2);
x = x(1);
% Derivatives of first basis function
Dxi{1} = [(2*x + 2*y + 1)*(2*x + 2*y - 2) + 2*(x + y - 1)^2, ...
          (x + y - 1)^2 + x*(2*x + 2*y - 2), ...
          y*(2*x + 2*y - 2); ...
          (2*x + 2*y + 1)*(2*x + 2*y - 2) + 2*(x + y - 1)^2, ...
          x*(2*x + 2*y - 2), ...
          (x + y - 1)^2 + y*(2*x + 2*y - 2)];
% Derivatives of second basis function
Dxi{2} = [-2*x*(2*x-3)-2*x^2, 2*x*(x-1)+x^2, 2*x*y; 0, 0, x^2];
% Derivatives of third basis function
Dxi{3} = [0, y^2, 0; -2*y*(2*y-3)-2*y^2, 2*x*y, 2*y*(y-1)+y^2];
% Derivatives of "bubble"
Dxi{4} = [-y*(x+y-1)-x*y; -x*(x+y-1)-x*y];
end