function Zred = reduceUnderApprox(Z,option,order)
% reduceUnderApprox - Reduces the order of a zonotope so that an
%                     under-approximation of the original set is obtained
%
% Syntax:  
%    Zred = reduceUnderApprox(Z,option)
%
% Inputs:
%    Z - zonotope object
%    option - used redction method (only 'sum' implemented so far)
%
% Outputs:
%    Zred - reduced zonotope
%
% Example: 
%    Z = zonotope(rand(2,10)-0.5*ones(2,10));
%
%    ZredSum = reduceUnderApprox(Z,'sum',3); 
%   
%    hold on
%    plot(Z,[1,2],'r');
%    plot(ZredSum,[1,2],'b');
%
% Other m-files required: none
% Subfunctions: see below
% MAT-files required: none
%
% See also: none

% Author:       Niklas Kochdumper
% Written:      19-November-2018
% Last update:  ---
% Last revision: ---

%------------- BEGIN CODE --------------

    if strcmp(option,'sum')
        Zred = reduceUnderApproxSum(Z,order);
    else
        error('Wrong value for input argument ''option''!');
    end
end

% Auxiliary Functions -----------------------------------------------------

function Zred = reduceUnderApproxSum(Z,order)

    % number of generators that stay unreduced
    n = size(Z.Z,1);
    N = floor(order*n) - 1;

    % select generators to reduce
    [c,G,Gred] = selectSmallestGenerators(Z,N);
    
    % under-approximate the generators that are reduced by one generator
    % corresponeding to the sum of generators
    g = sum(Gred,2);
    
    % construct the reduced zonotope object
    Zred = zonotope([c,G,g]);
end

function [c,G,Gred] = selectSmallestGenerators(Z,N)

    % obtain object properties
    c = Z.Z(:,1);
    G_ = Z.Z(:,2:end);

    % sort according to generator length
    temp = sum(G_.^2,1);
    [~,ind] = sort(temp,'descend');
    
    % split into reduce and unreduced generators
    G = G_(:,ind(1:N));
    Gred = G_(:,ind(N+1:end));
end


%------------- END OF CODE --------------