function K = rhct_der2(X,C,d)

% Calculate "stiffness" matrix \int D^2u' * C * D^2u dx on a reduced HCT
% element defined with the 3 points in the cols of X. C is given. NOTE:
% Setting C = [1 1 0; 1 1 0; 0 0 0] corresponds to the Laplacian since we
% use the convetion that D^2u = (u_xx, u_yy, u_xy)
%
% TAKEN FROM: https://www.degruyter.com/view/journals/cmam/12/4/article-p486.xml
%
% Do we trust this?

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
J = [E(:,3), -E(:,2)];
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
D{1} = zeros(3,3);
D{2} = zeros(3,3);
D{3} = zeros(3,3);
for k = 1:3
   [kp,km] = rotateindex(k);
   J = F(:,[kp,km]);
   H = [1, 0 0; [0;0], J.'];
   G = [J(2,2)^2, J(2,1)^2, -2*J(2,1)*J(2,2); ...
        J(1,2)^2, J(1,1)^2, -2*J(1,1)*J(1,2); ...
        -J(1,2)*J(2,2), -J(1,1)*J(2,1), J(1,1)*J(2,2)+J(2,1)*J(1,2)]./mu^2;
   for i = 1:length(gaussWeights)
       %Dxi = evalbasis(gaussPoints(i,:));
       D{k}  = Dxi{i}{1}*H*M{k};
       D{kp} = Dxi{i}{1}*H*M{kp} + Dxi{i}{2}*H + Dxi{i}{4}*b{k,kp}.';
       D{km} = Dxi{i}{1}*H*M{km} + Dxi{i}{3}*H + Dxi{i}{4}*b{k,km}.';
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
Dxi{1} = [-6+12*x(1)+12*x(2), -4+6*x(1)+4*x(2),    2*x(2); ...
          -6+12*x(1)+12*x(2),    2*x(1)       , -4+4*x(1)+6*x(2); ...
          -6+12*x(1)+12*x(2), -2+4*x(1)+2*x(2), -2+2*x(1)+4*x(2)];
Dxi{2} = [6-12*x(1), 6*x(1)-2, 2*x(2); 0, 0, 0; 0, 0, 2*x(1)];
Dxi{3} = [0, 0, 0; 6-12*x(2), 2*x(1), 6*x(2)-2; 0, 2*x(2), 0];
Dxi{4} = [-2*x(2); -2*x(1); 1-2*x(1)-2*x(2)];
end