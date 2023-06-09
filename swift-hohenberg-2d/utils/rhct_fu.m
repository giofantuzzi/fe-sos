function fu = rhct_fu(X,FUN,uDOF,d)

% Calculate "stiffness" matrix \int FUN(u) dx on a reduced HCT
% element defined with the 3 points in the cols of X.
%
% ADAPTED FROM: https://www.degruyter.com/view/journals/cmam/12/4/article-p486.xml

% ---------------------------------------------------------------------------- %
% Gauss quadrature
% ---------------------------------------------------------------------------- %
% cub2D = gaussData(d);
% gaussPoints = 0.5.*( cub2D(:,1:2) + 1 );
% gaussWeights = cub2D(:,3)./4;
[gaussPoints,gaussWeights] = PrecomputedGaussLeg2DTri(d);
persistent xi dold
if isempty(xi) || d~=dold
    dold = d;
    for i = 1:length(gaussWeights)
       xi{i} = evalbasis(gaussPoints(i,:));
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
fu = 0;
D = zeros(3,3);
for k = 1:3
   [kp,km] = rotateindex(k);
   H = [1, 0 0; [0;0], F(:,[kp,km]).'];
   for i = 1:length(gaussWeights)
       %xi = evalbasis(gaussPoints(:,i));
       tmp1 = xi{i}{1}*H;
       D(k,:)  = tmp1*M{k};
       D(kp,:) = tmp1*M{kp} + xi{i}{2}*H + xi{i}{4}*b{k,kp}.';
       D(km,:) = tmp1*M{km} + xi{i}{3}*H + xi{i}{4}*b{k,km}.';
       uinterp = [D(1,:), D(2,:), D(3,:)];
       uinterp(abs(uinterp)<1e-12) = 0;
       uinterp = uinterp*uDOF(:);
       fu = fu + (gaussWeights(i)*mu) .* FUN(uinterp);
   end
end

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
Dxi{1} = (1-x(1)-x(2))^2 .* [1+2*x(1)+2*x(2), x(1), x(2)] ;
Dxi{2} = x(1)^2 .* [3-2*x(1), x(1)-1, x(2)];
Dxi{3} = x(2)^2 .* [3-2*x(2), x(1), x(2)-1];
Dxi{4} = x(1).*x(2).*(1-x(1)-x(2));
end