function Zcubic = mixedCubicMultiplication(Z1,Z2,Z3,T,varargin)
% mixedCubicMultiplication - computes the set corresponding to the cubic 
%                            multiplication of three different zonotopes
%
% Description:
%   Calulates the following set:
%   { z = (x1' T x2) * x3 | x1 \in Z1, x2 \in Z2, x3 \in Z3 }
%
% Syntax:  
%    Zcubic = mixedCubicMultiplication(Z1,Z2,Z3,T)
%    Zcubic = mixedCubicMultiplication(Z1,Z2,Z3,T,ind)
%
% Inputs:
%    Z1 - zonotope object
%    Z2 - zonotope object
%    Z3 - zonotope object
%    T - third-order tensor
%    ind - cell-array containing the non-zero indizes of the tensor
%
% Outputs:
%    Zcubic - zonotope object representing the set of the cubic mapping
%
% Example: 
%    Z1 = zonotope(rand(2,4));
%    Z2 = zonotope(rand(2,5));
%    Z3 = zonotope(rand(2,3));
%
%    T{1,1} = rand(2);
%    T{1,2} = rand(2);
%    T{2,1} = rand(2);
%    T{2,2} = rand(2);
%
%    Zcub = mixedCubicMultiplication(Z1,Z2,Z3,T);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: cubicMultiplication, mixedMultiplication

% Author:       Niklas Kochdumper
% Written:      16-August-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% parse input arguments
if nargin >= 5
   ind = varargin{1};
else
   temp = 1:size(T,2);
   ind = repmat({temp},[size(T,1),1]);
end

% initialize variables
dim = length(ind);
N1 = size(Z1.Z,2);
N2 = size(Z2.Z,2);
N3 = size(Z3.Z,2);
Nq = N1*N2;

Zcub = zeros(dim,N1*N2*N3);

% loop over all system dimensions
for i = 1:length(ind)
   
   % loop over all quadratic matrices: \sum_k (x1' T_k x2) * x3_k 
   for k = 1:length(ind{i})

       % quadratic evaluation
       quadMat = Z1.Z' * T{i,ind{i}(k)} * Z2.Z;
       quadVec = reshape(quadMat,1,[]);   
       
       % multiply with Z3
       for j = 1:N3
          Zcub(i,(j-1)*Nq + 1 : j*Nq) = Zcub(i,(j-1)*Nq + 1 : j*Nq) + ...
                                        quadVec * Z3.Z(ind{i}(k),j);
       end
   end 
end

% construct the resulting zonotope
Zcubic = zonotope(Zcub);


%------------- END OF CODE --------------