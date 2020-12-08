function [cZquad] = mixedMultiplication(cZ1,cZ2,Q)
% mixedMultiplication - computes the quadratic map for two constrained 
%                       zontopes {x_1 * Q{i} * x_2 |x_1 \in cZ1, x_2 \in cZ2}
%
% Syntax:  
%    [cZquad] = mixedMultiplication(cZ1,cZ2,Q)
%
% Inputs:
%    cZ1 - conZonotope object
%    cZ2 - conZonotope object
%    Q - quadratic coefficients as a cell of matrices
%
% Outputs:
%    cZquad - resulting conZonotope object
%
% Example: 
%    cZ1 = conZonotope([0 1 0 1;0 1 2 -1],[-2 1 -1],2);
%    cZ2 = conZonotope([0 3 0 1;0 0 2 1],[1 0 1],1);
%
%    Q{1} = [1 2;-1 0];
%    Q{2} = [-2 -1;0 1];
%    cZquad = mixedMultiplication(cZ1,cZ2,Q);
%
%    plotFilled(cZ1,[1,2],'r');
%    figure
%    plotFilled(cZ2,[1,2],'b');
%    figure
%    plotFilled(cZquad,[1,2],'g');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: quadraticMultiplication, zonotope/mixedMultiplication

% Author:       Niklas Kochdumper
% Written:      13-August-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% rescale the constrained zonotopes to reduce the over-approximation of the
% quadratic map
if ~isempty(cZ1.A)
    cZ1 = rescale(cZ1,'iter');
end

if ~isempty(cZ2.A)
    cZ2 = rescale(cZ2,'iter');
end

% get generator matrices of the constrained zonotopes
Zmat1 = cZ1.Z;
Zmat2 = cZ2.Z;
dim = length(Q);

% initialize variables
N1 = size(Zmat1,2)-1;
N2 = size(Zmat2,2)-1;
c = zeros(dim,1);
G = zeros(dim,(N1+1)*(N2+1)-1);

% loop over all dimensions
for i = 1:dim
    
    % quadratic evaluation
    quadMat = Zmat1'*Q{i}*Zmat2;
    
    % center 
    c(i,1) = quadMat(1,1);
    
    % generators resulting from a multiplication with a zonotope center
    G(i,1:N2) = quadMat(1,2:end);
    G(i,N2+1:N1+N2) = quadMat(2:end,1);
    
    % transform quadratic matrix to vector
    quadVec = reshape(quadMat(2:end,2:end),1,[]);
    
    % remaining generators
    G(i, N2+N1+1:end) = quadVec;
end

% construct the new constraint matrix   
if isempty(cZ1.A) && isempty(cZ2.A)
    A = []; b = [];
else
    A = [cZ2.A, zeros(size(cZ2.A,1),(N1+1)*(N2+1)-N2-1);
         zeros(size(cZ1.A,1),N2), cZ1.A, zeros(size(cZ1.A,1),(N1+1)*(N2+1)-N2-N1-1)];

    b = [cZ2.b;cZ1.b];
end

% generate the resuling constrained zonotope
cZquad = conZonotope([c, G],A,b);


%------------- END OF CODE --------------