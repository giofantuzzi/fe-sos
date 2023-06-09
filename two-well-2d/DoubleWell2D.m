% 2D Minimization problem with double-well potential
%
% min_{H^1_0} \int_Omega 0.5*c * |grad(u)|^2 + (u-a)^2*(u-b)^2 dOmega
%
% Approximate global minimizer with sparse SOS optimization

% Clean up
clear
clc

% Parameters
hvals = 1/10;
epsilon_sq = 0.1;
omega=2;
chordalExtension = false;
shapes = 'square'; % Options: 'square', 'circle', 'ellipse'

% Loop
for i = 1:length(hvals)

    % Clear yalmip between iterations
    yalmip clear

    % Bound on DOF
    h = hvals(i);
    bound = 1/h;

    % Mesh and variables (no BCs)
    disp('Making mesh...')
    switch shapes
        case 'square'
            [nodes,elements,dirichlet]=MeshSquare(h);
        case 'circle'
            [nodes,elements,dirichlet]=MeshCircle(h);
        case 'ellipse'
            [nodes,elements,dirichlet]=MeshEllipse(h);
    end
    [ID,IEN,LM]=locator(nodes,elements,dirichlet);
    N=nnz(ID); %number of unknowns

    % Variables, no boundary ones
    % using msspoly variables from SPOT is faster
    alpha = msspoly('a',N);
    u = alpha(1)*ones(size(ID));
    u(logical(ID)) = bound*alpha;
    u(~logical(ID)) = dirichlet(:,2);
    u=u(:);

    % Cliques
    disp('Making cliques...')
    if chordalExtension
        cliques = [];
        opts.sparsemoment.order = 1;
    else
        cliques=LM(:,sum(spones(LM),1)==3);
        cliques=mat2cell(cliques',ones(size(cliques,2),1),3);
    end

    % Functional E, symbolic, using scaled variables to the unit box
    % Also, add pointwise bounds on the solution
    disp('Setting up objective & constraints (slow!)...')
    E = functional(u,nodes,elements,[-1,2,epsilon_sq]);
    [~,Epowers,Ecoeffs] = decomp(E);
    clear E alpha u

    % Recover in YALMIP
    disp('Recovering yalmip variables...')
    u = sdpvar(N,1);
    m = recovermonoms(Epowers,u);
    E = Ecoeffs*m;

    % Add inequalities
    g = 1- u.^2;

    % yalmip options
    opts = sdpsettings();
    opts.mosek.MSK_IPAR_INTPNT_ORDER_METHOD = 'MSK_ORDER_METHOD_APPMINLOC';

    % Solve sparse SOS problem
    disp('Calling solver!')
    [pstar,y,exponents,sol,model] = solvesparsemoment(u,E,[],g,omega,[],opts,cliques);

    % Plot solution
    clear usol
    usol(logical(ID)) = bound*extractlinearmoments(y, exponents);
    usol(~logical(ID))=dirichlet(:,2);
    T = trisurf(elements, nodes(:,1), nodes(:,2), usol);
    T.EdgeColor = 'interp';
    T.FaceColor = 'interp';
    view([0,90])
    caxis([-1.0 0])
    axis equal
    axis([0 1 0 1])
end