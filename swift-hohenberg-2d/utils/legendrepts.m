function [x, w] = legendrepts( n )
%LEGPTS    Legendre points and Gauss-Legendre quadrature weights.
%   LEGPTS(N) returns N Legendre points X in (-1,1).
%
%   [X, W] = LEGPTS(N) returns also a row vector W of weights for Gauss-Legendre
%   quadrature.

% Asymptotic formula (Tricomi) - only for positive x.
if ( mod(n,2) )
    s = 1;
else
    s = 0;
end
k = ((n+s)/2:-1:1).';
theta = pi*(4*k-1)/(4*n+2);
x = ( 1 - (n-1)/(8*n^3) - 1/(384*n^4)*(39-28./sin(theta).^2) ).*cos(theta);

% Initialise:
Pm2 = 1;
Pm1 = x;
PPm2 = 0;
PPm1 = 1;
dx = inf;
counter = 0;

% Loop until convergence:
while ( norm(dx, inf) > eps && counter < 10 )
    counter = counter + 1;
    for k = 1:n-1,
        P = ((2*k+1)*Pm1.*x-k*Pm2)/(k+1);
        Pm2 = Pm1;
        Pm1 = P;
        PP = ((2*k+1)*(Pm2+x.*PPm1)-k*PPm2)/(k+1);
        PPm2 = PPm1;
        PPm1 = PP;
    end
    % Newton step:
    dx = -P./PP;
    % Newton update:
    x = x + dx;
    % Reinitialise:
    Pm2 = 1;
    Pm1 = x;
    PPm2 = 0;
    PPm1 = 1;
end

% Once more for derivatives:
for k = 1:n-1,
    P = ( (2*k+1)*Pm1.*x - k*Pm2 ) / (k+1);
    Pm2 = Pm1;
    Pm1 = P;
    PP = ( (2*k+1)*(Pm2+x.*PPm1) - k*PPm2 ) / (k+1);
    PPm2 = PPm1;
    PPm1 = PP;
end

% Reflect for negative values:
x = [-x(end:-1:1+s) ; x];
ders = [PP(end:-1:1+s) ; PP];

% Quadrature weights:
w = 2./((1-x.^2).*ders.^2)';

end