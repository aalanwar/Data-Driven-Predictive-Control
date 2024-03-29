function obj = rA_inputTie(obj, options)
% rA_inputTie - simplified version of @linearSys > private > inputTie
%  (calculation without sets, only for obj.taylor.inputCorr)
%
% Syntax:  
%    options = rA_inputTie(obj, options)
%
% Inputs:
%    obj     - linearSys
%    options - options for the computation of reachable sets
%
% Outputs:
%    options - options for the computation of reachable sets
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none
%
% References: -
% 
% Author:       Mark Wetzlinger
% Written:      25-Aug-2019
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% load data from object structure
Apower = obj.taylor.powers;
E      = obj.taylor.error;
taylorTerms = options.taylorTerms;
r = options.timeStep;
dim = dimension(obj);

% initialize Asum
Asum_pos = zeros(dim);
Asum_neg = zeros(dim);

for i=2:(taylorTerms+1)
    % compute factor
    exp1 = -i/(i-1); exp2 = -1/(i-1);
    factor = (i^exp1-i^exp2)*options.factor(i); 
    
    % init Apos, Aneg
    Apos = zeros(dim);
    Aneg = zeros(dim);
    
    % obtain positive and negative parts
    pos_ind = Apower{i-1}>0;
    neg_ind = Apower{i-1}<0;
    
    Apos(pos_ind) = Apower{i-1}(pos_ind);
    Aneg(neg_ind) = Apower{i-1}(neg_ind);
    
    % compute powers; factor is always negative
    Asum_pos = Asum_pos + factor*Aneg; 
    Asum_neg = Asum_neg + factor*Apos;
end
% instantiate interval matrix
Asum = interval(Asum_neg,Asum_pos);

% compute error due to finite Taylor series according to internal document
% "Input Error Bounds in Reachability Analysis"
Einput = E*r;

%write to object structure
obj.taylor.inputF = Asum + Einput;
% ^... rewrite this equation when E is computed with the new method (?)

% compute vTrans 
vTrans = obj.B*options.uTrans;
% consider constant input
if ~isempty(obj.c)
    vTrans = vTrans + obj.c;
end
obj.taylor.inputCorr = obj.taylor.inputF * zonotope(vTrans);

end


%------------- END OF CODE --------------