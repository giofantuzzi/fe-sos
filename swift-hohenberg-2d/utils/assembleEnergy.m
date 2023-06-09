function E = assembleEnergy(x,y,TRI,u,r,b)

% Compute "energy"
% E(u) = \int |u_xx+u_yy|^2 - 2*(u_x^2 + u_y^2) + (1-r)*u^2 - b*u^3 + 0.5*u^4 dx

dquad = 12;
fun = @(u) (0.5.*u - b).*u.^3;
fprintf('\tAssembling matrices...\n')
[M, K, H] = make_fe_matrices(x,y,TRI,dquad);
Q = H - 2.*K + (1-r).*M;
E = u(1)*ones(size(u,2),1);
for i=1:size(u,2)
    fprintf('\tAssembling quadratic terms...\n')
    E(i) = E(i) + u(:,i).'*(Q*u(:,i));
    fprintf('\tAssembling nonlinear terms (super slow!)...\n')
    E(i) = E(i) + make_fu(fun,u(:,i),x,y,TRI,dquad);
end
E = E - u(1);

end
