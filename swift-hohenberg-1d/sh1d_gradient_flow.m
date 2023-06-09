% 1D Minimization problem for Swift-Hohenberg potential
%
% E(u) = \int |u_xx + u|^2 - r*u^2 - b*u^3 + 0.5*u^4 dx
%
% where the integral is over a domain [-L,L]. We fix parameters
%
%  r = -0.3
%  b = 1.2
%
% We enforce the boundary conditions
%
% u = u_x = 0 at x = L and x = -L
%
% and try to approximate minimizers using the gradient flow
%
% u_t = -u_xxxx - 2*u_xx - (1-r)*u + (3*b/2)*u^2 - u^3
%
% We reformulate this as
%
% u_t = v_xx - 2*u_xx - (1-r)*u + (3*b/2)*u^2 - u^3
% u_xx = -v
%
% and discretize this using a piecewise-linear FE discretization.

% Clean up
clear
clc
close all

% Swift-Hohenberg parameters
N = 2^8;            % number of elements
L = 32;             % domain half-period
r = 0.3;            
b = 1.2;
nIC = 5;          % number of random initial conditions to try
Dt = 1e-2;          % timestep
maxIter = 1e6;      % max iteration

% Mesh & random initial condition for gradient flow with 10 fourier modes
% Note: the IC has its boundary values removed
x = linspace(-L,L,N)';
nmodes = 10;
ampl = 4*rand(nmodes,nIC)-2;
u = sin(pi/L.*x(2:end-1)*(1:nmodes)) * ampl;
v = -(pi^2/L^2) .* sin(pi/L.*x*(1:nmodes)) * diag((1:nmodes).^2) * ampl;

% Mass matrix, stiffness matrix & handle to evaluate the nonlinearity
[M,K] = makeMassStiffness(x);
M0 = M(2:end-1,2:end-1);
K0 = K(2:end-1,2:end-1);
NNL = @(u) nonlinearTerms(x,u,r,b);

% LHS matrix
A = [(1+Dt*(1-r)).*M0-(2*Dt).*K0, Dt.*K(2:end-1,:); ...
    Dt.*K(:,2:end-1), -Dt.*M];

% Loop: A u_new = M0*u - Dt.*N(u)
% Cache cholesky factorization of M + c*Dt*K
[Rt,D,p] = ldl( A, 'lower', 'vector');
t = 0; iter = 0;
zz = sparse(N,nIC);
w = [u; v];
while iter < maxIter
    t = t + Dt;
    iter = iter + 1;
    RHS = [M0*w(1:N-2,:) + Dt.*NNL(w(1:N-2,:)); zz];
    w(p,:) = Rt.'\ ( D\ ( Rt\ RHS(p,:) ) );
    % Display and save progress every so many iterations
    if rem(iter,250)==0
        fprintf('iter = %i\n', iter);
        u = [zeros(1,nIC); w(1:N-2,:); zeros(1,nIC)];
        plot(x,u,'.-'); axis([-L,L,-1.5,2.5]);
        xlabel('x'), ylabel('u(x)')
        title(sprintf('t = %6.4f', t)), drawnow
    end
    if rem(iter,1000)==0; save('results/gf_latest.mat'); end
end
u = [zeros(1,nIC); w(1:N-2,:); zeros(1,nIC)];
v = w(N-1:end,:);

% Evaluate "energy" at the final timestep
E = sh_energy(x,u,v,r,b);
[E,i] = sort(E); u = u(:,i); v = v(:,i);
% Save data
save('results/gf_latest.mat')
