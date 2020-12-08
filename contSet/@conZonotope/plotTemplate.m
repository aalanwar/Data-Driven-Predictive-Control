function handle = plotTemplate(obj,varargin)
% plotTemplate - Plot an over-approximation of the 2D-projection of a 
%                constrained zonotope with a template polyhedron
%
% Syntax:  
%    handle = plotTemplate(obj)
%    handle = plotTemplate(obj,dim,numDir,plotOptions)
%
% Inputs:
%    obj - constrained zonotope object
%    dim - dimensions of the projection
%    numDir - number of directions for the template polyhedron
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
%    hold on
%    plot(cZono,[1,2],'b'); 
%    plotTemplate(cZono,[1,2],100,'g');
%    plotTemplate(cZono,[1,2],16,'r');
%    xlim([-2,0]);
%    ylim([-2,3]);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plot, plotFilled, plotFilledSplit, plotFilledTemplate

% Author:       Niklas Kochdumper
% Written:      31-October-2018
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

% default settings
dim = [1,2];
numDir = 16;
plotOptions = {'b'};

% parse input arguments
if nargin >= 2 && ~isempty(varargin{1})
   dim = varargin{1}; 
end
if nargin >= 3 && ~isempty(varargin{2})
   numDir = varargin{2}; 
end
if nargin >= 4 && ~isempty(varargin{3})
   plotOptions = varargin(3:end); 
end

% project the object to the 2D-subspace
obj = project(obj,dim);

% select directions for template polyhedron 
angles = linspace(0,360,numDir+1);
angles = angles(1:end-1);
angles = deg2rad(angles);

C = zeros(2,length(angles));

hold on
for i = 1:length(angles)
   C(:,i) = [cos(angles(i));sin(angles(i))];
end

% calculate the upper bounds along the directions
d = zeros(length(angles),1);

for i = 1:length(d)
    d(i) = boundDir(obj,C(:,i),'upper');
end

% construct template polyhedron
poly = mptPolytope(C',d);

% plot the template polyhedron
plot(poly,[1,2],plotOptions{:});
handle = gca;

%------------- END OF CODE --------------