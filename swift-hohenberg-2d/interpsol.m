function Up = interpsol(u,x,y,TRI,Xp,Yp)

% Interpolate solution from reduced HCT mesh onto given points with
% coordinates Xp and Yp

% Create matrix of points
[rows,cols] = size(Xp);
points = [Xp(:), Yp(:)].';
Up = zeros(size(points,2),1);

% Loop over elements
for i = 1:size(TRI,1)
    T = [x(TRI(i,:)), y(TRI(i,:))].';   % triangle
    mark = inTriangle(T,points);        % find points in this triangle
    Z = points(:,mark);                 % points to interpolate in this iter
    Up(mark) = rHCTinterp(u(:,TRI(i,:)),T,Z);
end
Up = reshape(Up,rows,cols);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function yn = inTriangle(T,X)
% Check if point with coordinates X is inside a triangle with coordinates T
% Can check m points at once (X is 2-by-m)
TOL = 1e-12;
A = T(:,2:3) - T(:,1);
s = A\(X-T(:,1));
yn = all(s>=-TOL,1) & all(s<=1+TOL,1) & all(sum(s)<=1+TOL,1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Up = rHCTinterp(u,X,Z)
% interpolate function on rHCT element given by the triangle T of the mesh

% Preliminary stuff
E = X(:,[3 1 2]) - X(:,[2 3 1]);        % Triangle edges
F = ( E(:,[2 3 1]) - E(:,[3 1 2]) )./3; % inter-element edges
R = [0 -1; 1 0];
N = R*E;
J = [E(:,2), E(:,3)];
mu = abs(det(J))/3;     % volume
a0 = sum(X,2)./3;

% vectors b{k,j}
b = cell(3,3);
for k = 1:3
    [kp,km] = rotateindex(k);
    normEk2 = dot(E(:,k),E(:,k));
    boo = 3*mu/normEk2;
    b{k,kp} = [6*dot(E(:,k), F(:,km))/normEk2; 2*F(:,km)+boo*N(:,k)];
    b{k,km} = [-6*dot(E(:,k), F(:,kp))/normEk2; 2*F(:,kp)+boo*N(:,k)];
end

% Matrices S, T and M
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

% Evaluate
Up = zeros(size(Z,2),1);
for k = 1:3
    [kp,km] = rotateindex(k);
    J = F(:,[kp,km]);
    s = J\(Z-a0);
    mark = find( all(s>=-1e-12,1) & all(s<=1+1e-12,1) & all(sum(s)<=1+1e-12,1) );
    H = [1, 0 0; [0;0], J.'];
    PSI = zeros(length(mark),9);
    for i = 1:length(mark)
        zloc = J\(Z(:,mark(i))-a0);
        Dxi = evalbasis(zloc);
        tmp = Dxi{1}*H;
        D{k}  = tmp*M{k};
        D{kp} = tmp*M{kp} + Dxi{2}*H + Dxi{4}*b{k,kp}.';
        D{km} = tmp*M{km} + Dxi{3}*H + Dxi{4}*b{k,km}.';
        PSI(i,:) = [D{1}, D{2}, D{3}];
    end
    Up(mark) = PSI*u(:);
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
Dxi{1} = (1-x(1)-x(2))^2 .* [1+2*x(1)+2*x(2), x(1), x(2)] ; % PHI0
Dxi{2} = x(1)^2 .* [3-2*x(1), x(1)-1, x(2)];    % PHI1
Dxi{3} = x(2)^2 .* [3-2*x(2), x(1), x(2)-1];    % PHI2
Dxi{4} = x(1).*x(2).*(1-x(1)-x(2));     % beta
end