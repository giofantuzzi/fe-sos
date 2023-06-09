function [M,K,D] = make_fe_matrices(x,y,TRI,dquad)

% Assemble FE matrices for the reduced HCT element. The list of DOF is a 
% 3-by-length(x) matrix (which we do not need to consider). We assemble:
% M: a "mass" matrix to discretize \int u^2 dx
% K: a "stiffness" matrix to discretize \int u_x^2 + u_y^2 dx
% D: a matrix to discretize the square laplacian, \int (u_xx + u_yy)^2 dx

% Initialize
nel = size(TRI,1);
iK = zeros(nel*81,1);
jK = zeros(nel*81,1);
vM = zeros(nel*81,1);
vK = zeros(nel*81,1);
vD = zeros(nel*81,1);

% Indexing of DOF
DOF = reshape(1:3*length(x),3,length(x));

% Loop
% C1 = [1 1 0; 1 1 0; 0 0 0]; % square of the Laplacian: (u_xx + u_yy)^2
C1 = [1 0 0; 0 1 0; 0 0 2]; % alternative formulation after integration by parts: u_xx^2 +2*u_xy^2 + u_yy^2
C2 = [1 0; 0 1]; % norm of the gradient: u_x^2 + u_y^2
idx = 1:81;
for i = 1:nel
    % Vertices of element & elemental matrix
    X = [x(TRI(i,:)), y(TRI(i,:))].';
    D = rhct_der2(X,C1,dquad); D(abs(D)<1e-12) = 0;
    K = rhct_der1(X,C2,dquad); K(abs(K)<1e-12) = 0;
    M = rhct_der0(X,dquad); M(abs(M)<1e-12) = 0;
    % Find row/col indices and store values
    localDOF = DOF(:,TRI(i,:)); 
    [c,r] = meshgrid(localDOF(:), localDOF(:));
    iK(idx) = r(:);
    jK(idx) = c(:);
    vD(idx) = D(:);
    vK(idx) = K(:);
    vM(idx) = M(:);
    % Update indices
    idx = idx + 81;
end

% Assemble
D = sparse(iK, jK, vD, numel(DOF), numel(DOF)); D(abs(D)<1e-12) = 0;
K = sparse(iK, jK, vK, numel(DOF), numel(DOF)); K(abs(K)<1e-12) = 0;
M = sparse(iK, jK, vM, numel(DOF), numel(DOF)); M(abs(M)<1e-12) = 0;
