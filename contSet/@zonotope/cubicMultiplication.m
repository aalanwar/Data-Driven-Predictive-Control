function Zcubic = cubicMultiplication(Z,T,varargin)
% cubicMultiplication - computes the set corresponding to the cubic 
%                       multiplication of a zonotope with a third-order
%                       tensor
%
% Description:
%   Calulates the following set:
%   { z = (x' T x) * x | x1 \in Z }
%
% Syntax:  
%    Zcubic = cubicMultiplication(Z,T)
%    Zcubic = cubicMultiplication(Z,T,ind)
%
% Inputs:
%    Z - zonotope object
%    T - third-order tensor
%    ind - cell-array containing the non-zero indizes of the tensor
%
% Outputs:
%    Zcubic - zonotope object representing the set of the cubic mapping
%
% Example: 
%    Z = zonotope(rand(2,4));
%
%    T{1,1} = rand(2);
%    T{1,2} = rand(2);
%    T{2,1} = rand(2);
%    T{2,2} = rand(2);
%
%    Zcub = cubicMultiplication(Z,T);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: quadraticMultiplication, mixedCubicMultiplication

% Author:       Niklas Kochdumper
% Written:      17-August-2018
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
N = size(Z.Z,2);
Zcub = zeros(dim,N^3);

% loop over all system dimensions
for i = 1:length(ind)
   
   listQuad = repmat({zeros(N,N)},[1,N]);
    
   % loop over all quadratic matrices: \sum_k (x' T_k x) * x_k 
   for k = 1:length(ind{i})

       % quadratic evaluation
       quadMat = Z.Z' * T{i,ind{i}(k)} * Z.Z;
       
       % add up all entries that correspond to identical factors
       temp = tril(quadMat,-1);
       quadMat = quadMat - temp;
       quadMat = quadMat + temp';
       
       % multiply with the zonotope generators of the corresponding
       % dimension
       for j = 1:N
          listQuad{j} = listQuad{j} + quadMat * Z.Z(ind{i}(k),j);
       end
   end 
   
   % add up all entries that belong to identical factors
   for k = 2:N
      
       % loop over all quadratic matrix rows whos factors already appear in
       % one of the previous quadratic matrices
       for j = 1:k-1
           
           % loop over all row entries
           for h = j:N
              if h <= k
                  listQuad{j}(h,k) = listQuad{j}(h,k) + listQuad{k}(j,h);
              else
                  listQuad{j}(k,h) = listQuad{j}(k,h) + listQuad{k}(j,h);
              end
           end
       end
   end
   
   % half the entries for purely quadratic factors
   temp = diag(listQuad{1});
   listQuad{1}(1,1) = listQuad{1}(1,1) + 0.5*sum(temp(2:end));
   
   for k = 2:N
      listQuad{1}(k,k) = 0.5 * listQuad{1}(k,k); 
   end
   
   % summerize all identical factors in one matrix
   counter = 1;
   
   for k = 1:N
   
       % loop over all matrix rows that contain unique factors
       for j = k:N
           m = N-j+1;         % number of elements in the row
           Zcub(i,counter : counter + m - 1) = listQuad{k}(j,j:end);
           counter = counter + m;
       end
   end
end

% concatenate the generator matrices
Zcub = Zcub(:,1:counter-1);

% construct the resulting zonotope
Zcubic = zonotope(Zcub);


%------------- END OF CODE --------------