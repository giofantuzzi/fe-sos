function D = makeCliques(u,TRI,isbnd)

% Assemble cliques (work by elements!)
nel = size(TRI,1);
D = cell(nel,1);
numnodes = sum(~isbnd);
DOF = zeros(size(u));
DOF(:,~isbnd) = reshape(1:3*numnodes,3,numnodes);
for i = 1:nel
    % Do we hit the boundary?
    if ~any(isbnd(TRI(i,:)))
            D{i} = DOF(:,TRI(i,:));
            D{i} = sort(D{i}(:)).';
    end
end

% Remove empty ones
keep = cellfun(@(X)numel(X)>0, D, 'UniformOutput', 1);
D = D(keep);