function B = assembleHermiteMatrix(MESH,DER,CLEANTOL)

%initialise
num_elements = length(MESH)-1;
Q = hermiteIntegrationMatrix(0,1,DER); % Standard matrix, for sizing
[nQ,mQ] = size(Q); % nQ should be even, mQ even or 1
idx1 = 1:nQ;
idx2 = 1:mQ;
num_terms = nQ*mQ;

% shifts
shift_idx1 = nQ/2;
shift_idx2 = floor(mQ/2); % either mQ/2, or 0 if mQ==1

% Containers
iB = zeros(num_terms*num_elements,1);
jB = zeros(num_terms*num_elements,1);
vB = zeros(num_terms*num_elements,1);

% Loop
shift = 0;
for i=1:num_elements
    Q = hermiteIntegrationMatrix(MESH(i),MESH(i+1),DER);
    pos = shift + (1:num_terms);
    [c,r] = meshgrid(idx2,idx1);
    iB(pos) = r(:);
    jB(pos) = c(:);
    vB(pos) = Q(:);
    shift = shift + num_terms;
    idx1 = idx1 + shift_idx1;
    idx2 = idx2 + shift_idx2;
end

% Size of B
nNodes = num_elements+1;
nB = shift_idx1*nNodes;
mB = 1;
if mQ>1
    mB = shift_idx2*nNodes;
end

% Assemble and clean
B = sparse(iB, jB, vB, nB, mB);
[iB, jB, vB] = find(B);
keep = abs(vB)>=CLEANTOL;
B = sparse(iB(keep), jB(keep), vB(keep), nB, mB);

% End function
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% UTILITY FUNCTION
function Q = hermiteIntegrationMatrix(a,b,DER)

% Order correctly
if a>b
    c = b;
    b = a;
    a = c;
elseif a==b
    error('You must specify a < b!')
end

if length(DER)==1
    % LINEAR TERMS
    if DER==0
        % int_a^b eta dx
        Q = [ (b-a)/2, (a - b)^2/12, (b-a)/2, -(a - b)^2/12].';
    elseif DER==1
        % int_a^b Deta dx
        Q = [ -1, 0, 1, 0].';
    elseif DER==2
        % int_a^b D2eta dx
        Q = [ 0, -1, 0, 1].';
    end
    
elseif length(DER)==2
    % QUADRATIC TERMS
    if DER(1)==0 && DER(2)==0
        % int_a^b eta*eta.' dx
        Q = [ 13*(b-a)/35, (11*(a - b)^2)/210,   9*(b-a)/70, -(13*(a - b)^2)/420;...
            (11*(a - b)^2)/210,     -(a - b)^3/105,    (13*(a - b)^2)/420,       (a - b)^3/140;...
            9*(b-a)/70, (13*(a - b)^2)/420, 13*(b-a)/35, -(11*(a - b)^2)/210;...
            -(13*(a - b)^2)/420,      (a - b)^3/140,   -(11*(a - b)^2)/210,      -(a - b)^3/105];
    elseif DER(1)==1 && DER(2)==1
        % int_a^b Deta*Deta.' dx
        Q = [ -6/(5*(a - b)), 1/10,  6/(5*(a - b)), 1/10;...
               1/10, (2*b)/15 - (2*a)/15, -1/10, (a-b)/30;...
               6/(5*(a - b)), -1/10, -6/(5*(a - b)), -1/10;...
               1/10, (a-b)/30, -1/10, (2*b)/15 - (2*a)/15];
    elseif DER(1)==2 && DER(2)==2
        % int_a^b D2eta*D2eta.' dx
        Q = [ -12/(a - b)^3,  6/(a - b)^2,  12/(a - b)^3,  6/(a - b)^2;...
             6/(a - b)^2,   -4/(a - b),  -6/(a - b)^2,   -2/(a - b);...
             12/(a - b)^3, -6/(a - b)^2, -12/(a - b)^3, -6/(a - b)^2;...
             6/(a - b)^2,   -2/(a - b),  -6/(a - b)^2,   -4/(a - b)];
    end

elseif length(DER)==3
   % CUBIC TERMS
   if DER(1)==0 && DER(2)==0 && DER(3)==1
   Q = [-1/3,     (5*b)/42 - (5*a)/42,                     1/3, (17*a)/210 - (17*b)/210;...
        (5*a)/84 - (5*b)/84,           (a - b)^2/168,     (5*b)/84 - (5*a)/84,     -(11*(a - b)^2)/840;...
        -1/6,   (2*a)/105 - (2*b)/105,                     1/6,   (2*a)/105 - (2*b)/105;...
        (17*b)/420 - (17*a)/420,           (a - b)^2/280, (17*a)/420 - (17*b)/420,           (a - b)^2/168;...
        (5*a)/84 - (5*b)/84,           (a - b)^2/168,     (5*b)/84 - (5*a)/84,     -(11*(a - b)^2)/840;...
        -(a - b)^2/84,                       0,            (a - b)^2/84,           (a - b)^3/420;...
        (17*a)/420 - (17*b)/420,          -(a - b)^2/168, (17*b)/420 - (17*a)/420,          -(a - b)^2/280;...
        (a - b)^2/105,          -(a - b)^3/840,          -(a - b)^2/105,          -(a - b)^3/840;...
        -1/6,   (2*a)/105 - (2*b)/105,                     1/6,   (2*a)/105 - (2*b)/105;...
        (17*a)/420 - (17*b)/420,          -(a - b)^2/168, (17*b)/420 - (17*a)/420,          -(a - b)^2/280;...
        -1/3, (17*a)/210 - (17*b)/210,                     1/3,     (5*b)/42 - (5*a)/42;...
        (5*b)/84 - (5*a)/84,      (11*(a - b)^2)/840,     (5*a)/84 - (5*b)/84,          -(a - b)^2/168;...
        (17*b)/420 - (17*a)/420,           (a - b)^2/280, (17*a)/420 - (17*b)/420,           (a - b)^2/168;...
        (a - b)^2/105,          -(a - b)^3/840,          -(a - b)^2/105,          -(a - b)^3/840;...
        (5*b)/84 - (5*a)/84,      (11*(a - b)^2)/840,     (5*a)/84 - (5*b)/84,          -(a - b)^2/168;...
        -(a - b)^2/84,           (a - b)^3/420,            (a - b)^2/84,                       0];
   end
   
end

end

%%%%% SYMBOLIC COMPUTATIONS TO FIND MORE %%%%%
% clear all
% clc
% syms pa pb Dpa Dpb
% syms a b x
% Q = [1 a a^2 a^3; 1 b b^2 b^3; 0 1 2*a 3*a^2; 0 1 2*b 3*b^2];
% Q = simplify(Q\[pa;pb;Dpa;Dpb]);
% p = simplify(Q.'*[1; x; x^2; x^3]);
% [xi,mmm] = coeffs(p,[pa;Dpa;pb;Dpb]);
% xi =  simplify(xi.');
% Dxi = simplify(diff(xi,x));
% D2xi = simplify(diff(Dxi,x));
% Q00 = simplify( int(xi*xi.', x, a, b) )
% Q11 = simplify( int(Dxi*Dxi.', x, a, b) )
% Q22 = simplify( int(D2xi*D2xi.', x, a, b) )
% eta = [(b-x)/(b-a); (x-a)/(b-a)];
% Deta = diff(eta,x);
% M = xi*xi.'; M = M(:);
% N = simplify( int(M*Deta.', x, a, b) )
