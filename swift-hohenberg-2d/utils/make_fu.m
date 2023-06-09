function fu = make_fu(FUN,uDOF,x,y,TRI,dquad)

% Quadrature of FUN(u) over mesh of HCT elements. The input FUN is a
% function handle that can work on symbolic variables, while uDOF is a
% vector with 3*length(x) degrees of freedom.

% Initialize
uDOF = reshape(uDOF,3,length(x));

% First element (initalize of the same kind as uDOF, so it works with
% yalmip variables without throwing NaNs)
X = [x(TRI(1,:)), y(TRI(1,:))].';
fu = rhct_fu(X, FUN, uDOF(:,TRI(1,:)), dquad);

% Loop over remaining element
if isa(uDOF,'sdpvar')
    for i = 2:size(TRI,1)
        % Vertices of element & elemental matrix
        X = [x(TRI(i,:)), y(TRI(i,:))].';
        K = rhct_fu(X, FUN, uDOF(:,TRI(i,:)), dquad);
        fu = fu + K;
    end
else
    parfor i = 2:size(TRI,1)
        % Vertices of element & elemental matrix
        X = [x(TRI(i,:)), y(TRI(i,:))].';
        K = rhct_fu(X, FUN, uDOF(:,TRI(i,:)), dquad);
        fu = fu + K;
    end
end