function [M, K, H] = makeMassStiffness(x,type)

if nargin<2; type = []; end

if isempty(type) ||  strcmpi(type,'lagrange')
    % Piecewise linear Lagrange interpolation
    % Mass and stiffness matrix on standard element [-1,1]
    M0 = [2 1 1 2]'./3;
    K0 = [1 -1 -1 1]./2;
    
    % Initialize
    N = length(x)-1;
    i = zeros(4*N,1);
    j = zeros(4*N,1);
    vM = zeros(4*N,1);
    vK = repmat(4*N,1);
    
    % Loop over elements to construct sparse matrix
    idx = 1:2;
    pos = 1:4;
    for k = 1:N
        % row/col indices
        [c,r] = meshgrid(idx,idx);
        i(pos) = r(:);
        j(pos) = c(:);
        % values (take mesh scaling into accound)
        s = 2./( x(k+1) - x(k) );
        vM(pos) = M0./s;
        vK(pos) = K0.*s;
        % Update pos and idx
        idx = idx + 1;
        pos = pos + 4;
    end
    
    % Create matrices
    M = sparse(i,j,vM,N+1,N+1);
    K = sparse(i,j,vK,N+1,N+1);
    H = [];
    
elseif strcmpi(type,'hermite')
    % Cubic hermite interpolation
    M = assembleHermiteMatrix(x,[0 0],1e-12);
    K = assembleHermiteMatrix(x,[1 1],1e-12);
    H = assembleHermiteMatrix(x,[2 2],1e-12);
else
    error('Type of element not recognized')
end
end
