function [V,H,Hlast,happyBreakdown] = arnoldi(A,vInit,redDim)
% arnoldi - computes the Arnoldi iteration for building a Krylov subspace
%
% Syntax:  
%    H = arnoldi(obj,options)
%
% Inputs:
%    obj - linearSys object
%    vInit - initial value of the vector that is multiplied
%    options - reachability options
%
% Outputs:
%    H - transformation matrix
%
% Example: 
%    -
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Matthias Althoff
% Written:      15-November-2016 
% Last update:  22-December-2016
%               06-November-2018
% Last revision:---

%------------- BEGIN CODE --------------

%preallocate H
H = zeros(redDim,redDim);

% initialize 
V(:,1) = vInit/norm(vInit);
happyBreakdown = 0;
    
% compute elements of transformation matrix
for j = 1:redDim
    % update v
    w = A*V(:,j);
    
    % generate column vector of H
    for i = 1:j
        H(i,j) = w'*V(:,i);
        w = w - H(i,j)*V(:,i);
    end
    
    % last element and update of q{k}
    H(j+1,j) = norm(w);
    % happy-breakdown?
    if H(j+1,j) <= eps %*norm_A
        happyBreakdown = 1;
        break
    end
    V(:,j+1) = w/H(j+1,j);
end    

% save H(j+1,j)
Hlast = H(j+1,j);

% remove last column of V
if ~happyBreakdown % no happy breakdown
    V(:,j+1) = [];
    % remove last row of H
    H(j+1,:) = [];
else
    % reduce H due to happy breakdown
    Htmp = H;
    H = Htmp(1:j,1:j);
end


%------------- END OF CODE --------------