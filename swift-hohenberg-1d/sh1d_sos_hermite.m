% 1D Minimization problem for Swift-Hohenberg potential
%
% E(u) = \int |u_xx + u|^2 - r*u^2 - b*u^3 + 0.5*u^4 dx
%
% where the integral is over a domain [-L,L]. We fix parameters
%
%  r = 0.3
%  b = 1.2
%
% We enforce the boundary conditions
%
% u = u_x = 0 at x = L and x = -L
%
% and try to approximate minimizers using sparse moment-SOS relaxations with a
% mixed formulation. That is, we introduce a "slack" variable v subject to
%
% <w_x, u_x> = <w,v> for all w \in H^1(-L,L)
%
% which implies that u is in H^2_0(-L,L). This formulation results in a POP
% with equality constraints

% Clean up
clear
clc
yalmip clear

% Swift-Hohenberg parameters
L = 32;
r = 0.3;
b = 1.2;

% Discretization parameters (2^pow elements)
pow = 5;

% Moment-SOS relaxation parameters
omega = 4;

% Loop
for i = 1:length(pow)
    
    yalmip clear
    
    % Mesh and variables (no BCs) -- use SPOT for symbolic assembly, much
    % faster than YALMIP
    N = 2^pow(i)+1;
    x = linspace(-L,L,N);
    u = msspoly('u',2*N-4);
    uwithbc = [0; 0; u; 0; 0];
    
    % Bound
    h = L*2/(N-1);
    bound = 2/h;
    
    % Set up the POP using scaled variables to the unit box
    % Take into account the boundary conditions on u
    disp('Setting up objective & constraints...')
    E = shEnergy(x,bound*uwithbc,[],r,b,'hermite');
    [~,Epowers,Ecoeffs]=decomp(E);
    
    % Recover in YALMIP
    disp('Recovering yalmip variables...')
    u = sdpvar(2*N-4,1);
    m = recovermonoms(Epowers,u);
    E = Ecoeffs*m;
    
    % Inequality constraints (bounds on DOFs)
    g = 1-u.^2;
    
    % Solve sparse SOS problem (chordal by construction so can safely use
    % automatic detection of cliques) Then, attempt to extract solution by reading
    % the degree-1 moments
    disp('Calling solver!')
    opts = sdpsettings();
    opts.sparsemoment.mergeCliques = 1;
    [pstar,y,exponents,sol,model] = solvesparsemoment(u,E,[],g,omega,[],opts);
    
    % Attempt to extract solution
    try
        usol = extractminimizers(sol.momentMatrices,sol.gramMonomials);
        usol = matchsol(sol.cliques.Set,usol);
    catch
        fprintf('Minimizer extraction failed! Resorting to linear moments...\n')
        usol = extractlinearmoments(y, exponents);
    end
    
    % Plot approximate solution using 50 points for each element
    for mm = 1:size(usol,2)
        vv = bound*[0; 0; usol(:,mm); 0; 0];
        plotHermite(x,vv,50,'.-')
        hold on 
        plotHermite(x,vv,2,'o')
        xlim([-L,L])
        if mm > 1; hold on; end
        drawnow
    end
    hold off
end