function [cZquad] = quadraticMultiplication(cZ, Q)
% quadraticMultiplication - computes the quadratic map 
%                           {x^T Q{i} x | x \in cZ}
%
% Syntax:  
%    [cZquad] = quadraticMultiplication(cZ, Q)
%
% Inputs:
%    cZ - conZonotope object
%    Q - quadratic coefficients as a cell of matrices
%
% Outputs:
%    cZquad - conZonotope object
%
% Example: 
%    Z = [0 1 0 1;0 1 2 -1];
%    A = [-2 1 -1];
%    b = 2;
%    cZ = conZonotope(Z,A,b);
%
%    Q{1} = [1 2;-1 0];
%    Q{2} = [-2 -1;0 1];
%    cZquad = quadraticMultiplication(cZ,Q);
%
%    plotFilled(cZ,[1,2],'r');
%    figure
%    plotFilled(cZquad,[1,2],'b');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mixedMultiplication, zonotope/quadraticMultiplication

% Author:       Niklas Kochdumper
% Written:      13-August-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% rescale the constrained zonotope to reduce the over-approximation of the
% quadratic map
if ~isempty(cZ.A)
    cZ = rescale(cZ,'iter');
end

% get generator matrix of constrained zonotope
Zmat = cZ.Z;
dim = length(Q);
gens = length(Zmat(1,:)) - 1;

% initialize variables
N = (gens+1) * 0.5*(gens+2)-1;
c = zeros(dim,1);
G = zeros(dim,N);

% loop over all dimensions
for i = 1:dim
    
    % quadratic evaluation
    quadMat = Zmat'*Q{i}*Zmat;
    
    % center
    c(i,1) = quadMat(1,1) + 0.5*sum(diag(quadMat(2:end,2:end)));
    
    % generators from diagonal elements
    ind = 1:gens;
    G(i, ind) = 0.5*diag(quadMat(ind+1,ind+1));
    
    % generators from other elements
    counter = 0;
    for j = 0:gens
        kInd = (j+1):gens;
        G(i, gens + counter + kInd - j) = quadMat(j+1, kInd+1) + quadMat(kInd+1, j+1)';
        counter = counter + length(kInd);
    end
end

% construct the constrained matrices for the constrained zonotope that
% represents the resuling set for the quadratic map
N_ = N - 2*gens;
A = [zeros(size(cZ.A,1),gens), cZ.A, zeros(size(cZ.A,1),N_)];

% generate the resuling conZonotope object
cZquad = conZonotope([c, G],A,cZ.b);


%------------- END OF CODE --------------