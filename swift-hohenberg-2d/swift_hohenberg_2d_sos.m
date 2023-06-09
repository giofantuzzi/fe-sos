% Minimization problem for 2D Swift-Hohenberg potential
%
% E(u) = \int |u_xx + u_yy + u|^2 - r*u^2 - b*u^3 + 0.5*u^4 dx
%
% where the integral is over a domain [-2L,2L] x [-L,L]. We fix parameters
%
%  r = 0.3
%  b = 1.2
%
% We enforce boundary conditions of H^2_0 type and try to approximate minimizers
% using sparse moment-SOS relaxations. We use a reduced 9-DOF HCT element!

% Clean up
clear

% Swift-Hohenberg parameters
L = 6;
r = 0.3;
b = 1.2;

% Discretization parameters
% 2*N elements in x, N elements in y
Nvals = 10; 

% Moment-SOS relaxation parameters
% (TOO EXPENSIVE FOR omega > 2)
omega = 2;

% Loop
for i = 1:length(Nvals)
        
    % Mesh
    N = Nvals(i);
    x = linspace(-2*L, 2*L, 2*N+1);
    y = linspace(-L, L, N+1);
    [X,Y] = meshgrid(x,y);
    x = X(:);
    y = Y(:);
    TRI = delaunay(x,y);
    isbnd = (abs(x-2*L)<1e-12) | (abs(x+2*L)<1e-12);
    isbnd = isbnd | (abs(y-L)<1e-12) | (abs(y+L)<1e-12);
    
    % Bounds
    h = sqrt(2)*2*L/N;
    bound = 4; 
    
    % Load problem data from ./input-data/ if it exists
    FNAME = sprintf('input-data/L%i/poly_data_N%i.mat',L,N);
    if exist(FNAME,'file')
        load(FNAME)
    else
        % Add utils to set up problem with SPOT
        % This is slow!
        addpath(genpath('./utils'))

        % Spot variables
        u = reshape(msspoly('u',3*sum(~isbnd)), 3, sum(~isbnd));
        uwithbc = u(1)*ones(3,length(x));
        uwithbc(:,find(~isbnd)) = u;
        uwithbc(:,find(isbnd)) = 0;
        u = u(:);
        
        % "Cliques" (they do NOT satisfy the RIP)
        cliques = makeCliques(uwithbc,TRI,isbnd);
        
        % Set up the POP using unscaled variable -- will be scaled later
        % when recovering the problem for optimization.
        % Take into account the boundary conditions on u
        disp('Setting up objective & constraints (slow!)...')
        E = assembleEnergy(x,y,TRI,uwithbc(:),r,b);
        [~,Epowers,Ecoeffs]=decomp(E);
        keep = abs(Ecoeffs)>1e-10;
        Ecoeffs = Ecoeffs(keep);
        Epowers = Epowers(keep,:);
        
        % Save data and clear useless variables
        save(FNAME,'Epowers','Ecoeffs','cliques')
        clear E u uwithbc
        
        % Clean path
        rmpath(genpath('./utils'))
    end
    
    % Recover the problem, scaling the variable according to the specified
    % bound so the actual optimization variable in the POP lies in a unit
    % box or ball (good for numerical conditioning)
    disp('Recovering problem...')
    nvars = 3*sum(~isbnd);
    SF = prod(bound.^Epowers, 2); % scale each DOF by "bound" in POP
    E.pows = Epowers; % powers don't change
    E.coef = Ecoeffs(:) .* SF(:); % but pick up a scaling factor for each monomial!
    scale = mean(abs(E.coef)); % rescale
    E.coef = E.coef ./ scale;
    clear Epowers Ecoeffs SF % These are no longer needed
    
    % Add inequalities
    % l2 norm of DOF in an element is bounded
    g = [];
    for k = 1:length(cliques)
        nvar_clique = length(cliques{k});
        g(k).pows = sparse(2:nvar_clique+1, cliques{k}, 2, nvar_clique+1, nvars);
        g(k).coef = [1; -ones(nvar_clique,1)];
    end
    
    % Set options for MOSEK (useful for HPC cluster)
    opts = sdpsettings();
    opts.mosek.MSK_IPAR_INTPNT_ORDER_METHOD = 'MSK_ORDER_METHOD_APPMINLOC';
    opts.mosek.MSK_IPAR_NUM_THREADS = 6;
    
    % Solve sparse SOS problem (chordal by construction so can safely use
    % automatic detection of cliques)
    disp('Calling solver!')
    [pstar,mom,exponents,sol] = solvesparsemoment_fast(nvars,E,[],g,omega,1,opts,cliques);
    pstar = pstar * scale;
    mosekTime = sol.solverTime;

    % Attempt to extract solution
    try
        xsol = extractminimizers(sol.momentMatrices,sol.gramMonomials);
        xsol = matchsol(sol.cliques.Set,xsol);
    catch
        fprintf('Minimizer extraction failed! Resorting to linear moments...\n')
        xsol = extractlinearmoments(mom, exponents);
    end
    if size(xsol,2)>1
        fprintf('Too many optimizers -- something wrong? Resorting to linear moments...\n')
        xsol = extractlinearmoments(mom, exponents);
    end
    
    % Interpolate for plotting
    Np = 2^7;
    [Xp,Yp] = meshgrid(linspace(-2*L, 2*L, 2*Np+1),linspace(-L, L, Np+1));
    for k = 1:min(size(xsol,2),3)
        usol = zeros(3,length(x));
        usol(:,~isbnd) = reshape(bound.*xsol(:,k), 3, sum(~isbnd));
        Up = interpsol(usol,x,y,TRI,Xp,Yp);
    end
    plotsolution(Xp,Yp,Up,TRI,x,y,0)
end