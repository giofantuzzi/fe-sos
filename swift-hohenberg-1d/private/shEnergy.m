function E = shEnergy(x,u,v,r,b,type)

% Compute the swift-hohenberg energy functional:
% E(u) = \int |u-v|^2 - r*u^2 - b*u^3 + 0.5*u^4 dx

if nargin<6; type = []; end

% Nonquadratic terms in nested function at the end of the script
fun = @(u) nnl(u,b);

% Initialize to ensure everything works -- this will be overwritten anyway
E = sum(u,1);

% Compute
if isempty(type) || strcmpi(type,'lagrange')
    [M, ~] = makeMassStiffness(x,'lagrange');
    umv = u - v;
    for i=1:size(u,2)
        E(i) = umv(:,i).'*(M*umv(:,i));
        E(i) = E(i) - r.*( u(:,i).'*(M*u(:,i)) );
        E(i) = E(i) + gaussQuadrature(fun,x,u(:,i),type);
    end
    
elseif strcmpi(type,'hermite')
    [M, K, H] = makeMassStiffness(x,'hermite');
    Q = H - 2.*K + (1-r).*M;
    for i=1:size(u,2)
        E(i) = u(:,i).'*(Q*u(:,i));
        E(i) = E(i) + gaussQuadrature(fun,x,u(:,i),type);
    end
else
    error('Type of element not recognized')
    
end

end

%%%%%%%%%%%%%%%%%%%%%
function v = nnl(u,b)
u2 = u.^2;
v = u2.*( 0.5*u2 - b*u );
end