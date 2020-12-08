function handle = plotFilledSplit(obj,varargin)
% plotFilledSplit - Plot an approximation of the 2D-projection of a 
%                   constrained zonotope object by applying recursive
%                   splits to the set
%
% Syntax:  
%    handle = plotFilledSplit(obj)
%    handle = plotFilledSplit(obj,dim,splits,plotOptions)
%
% Inputs:
%    obj - constrained zonotope object
%    dim - dimensions of the projection
%    splits - number of recursive set splits to get better approximation
%    plotOptions - plot settings specified as name-value pairs
%
% Outputs:
%    handle - object handle for the resulting graphic object
%
% Example: 
%    Z = [0 1 0 1;0 1 2 -1];
%    A = [-2 1 -1];
%    b = 2;
%    cZono = conZonotope(Z,A,b);
%
%    figure
%    plotFilledSplit(cZono,[1,2],4,'b','EdgeColor','none');
%
%    figure
%    plotFilledSplit(cZono,[1,2],8,[0,0.5,0],'EdgeColor','none');
%
%    figure
%    plotFilled(cZono,[1,2],'r','EdgeColor','none');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plot, plotFilled, plotFilledTemplate

% Author:       Niklas Kochdumper
% Written:      31-October-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% default settings
dim = [1,2];
splits = 4;
plotOptions = {'b'};

% parse input arguments
if nargin >= 2 && ~isempty(varargin{1})
   dim = varargin{1}; 
end
if nargin >= 3 && ~isempty(varargin{2})
    splits = varargin{2};
end
if nargin >= 4 && ~isempty(varargin{3})
   plotOptions = varargin(3:end); 
end

% project the object to the 2D-subspace
obj = project(obj,dim);

% recursively split the constrained zonotope
list = {obj};

for i = 1:splits
   
    listTemp = cell(length(list)*2,1);
    counter = 1;
    
    % loop over all sets at the current recursion level
    for j = 1:length(list)
        
       % calculate radius of the interval over-approximation as a heuristic
       % indicating which dimension should be best splitted
       inter = interval(list{j});
       r = rad(inter);
        
       % split the set
       if r(1) > r(2)
           temp = split(list{j},1);
       else
           temp = split(list{j},2);
       end
       
       % update variables
       listTemp{counter} = temp{1};
       listTemp{counter+1} = temp{2};
       counter = counter + 2;
    end
    
    list = listTemp;
end

% over-approximate the splitted sets with intervals
hold on

for i = 1:length(list)
   plotFilled(interval(list{i}),[1,2],plotOptions{:}); 
end

handle = gca;

%------------- END OF CODE --------------