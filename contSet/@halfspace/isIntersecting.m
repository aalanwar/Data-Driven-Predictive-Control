function res = isIntersecting(obj1,obj2)
% isIntersecting - determines if hyperplane obj1 intersects obj2
%
% Syntax:  
%    res = isIntersecting(obj1,obj2)
%
% Inputs:
%    obj1 - halfspace object
%    obj2 - interval or zonotope object
%
% Outputs:
%    result - 1/0 if set is intersecting, or not
%
% Example: 
%    ---
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: ---

% Author:       Niklas Kochdumper
% Written:      16-May-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

res = 1;

% convert interval to zontope
if isa(obj2,'interval')
   obj2 = zonotope(obj2); 
end

% different set representations
if isa(obj2,'conZonotope')
    
    % compute bounds along the hyperplane normal vector
    bu = boundDir(obj2,obj1.c,'upper');
    bl = boundDir(obj2,obj1.c,'lower');
    
    % check for intersection
    if obj1.d < bl|| obj1.d > bu
       res = 0; 
    end
    
elseif isa(obj2,'zonotope')

    % project zonotope center and generators onto hyperplane normal vector
    c_ = obj1.c' * obj2.Z(:,1);
    G_ = obj1.c' * obj2.Z(:,2:end);

    sumG = sum(abs(G_));

    % check for intersection
    if obj1.d < (c_ - sumG) || obj1.d > (c_ + sumG)
       res = 0; 
    end

else
   error('Operation not implemented!'); 
end


%------------- END OF CODE --------------