function NNL = nonlinearTerms(x,u,r,b)

% Shape functions, quadrature points and weights
% Use 5-point quadrature
gnodes = [-1; -sqrt(3/7); 0; sqrt(3/7); 1];
gweights = [1/10 49/90 32/45 49/90 1/10];
psi = 0.5.*[1-gnodes, gnodes+1]; % basis functions at the gauss nodes

% augment u with BCs
u = [zeros(1, size(u,2)); u; zeros(1,size(u,2))];

% Initialize
N = length(x)-1;
NNL = zeros(size(u));

% Loop over elements to construct sparse matrix
idx = 1:2;
bscaled = 1.5*b;
fun = @(u) (bscaled - u).*(u.*u);
for k = 1:N
    % quadrature of (1-r)*u^2 - b*u^3 + 0.5*u^4 on standard element
    % uloc is a column vector with values of u at the quadrature points
    s = 2./( x(k+1) - x(k) );
    uloc = psi*u(idx,:);
    FUNVAL = fun(uloc);
    tmp1 = FUNVAL.*psi(:,1);
    tmp2 = FUNVAL.*psi(:,2);
    NNL(idx,:) = NNL(idx,:) + [gweights*tmp1; gweights*tmp2]/s;
    % Update idx
    idx = idx + 1;
end

% Remove boundary points
NNL = NNL(2:end-1,:);
end