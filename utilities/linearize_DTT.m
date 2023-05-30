function [f0,A_lin,B_lin] = linearize_DTT(sysDisc,xstar,ustar)
% linearize - linearizes the nonlinearSysDT object
%
% Syntax:  
%    [obj,A_lin,U] = linearize(obj,R,options)
%
% Inputs:
%    obj - nonlinearSysDT system object
%    R - initial reachable set
%    options - options struct
%
% Outputs:
%    obj - nonlinearSysDT system object with additional properties
%    A_lin - system matrix of the linearized system
%    U - reachable set due to the inputs
%
% Example: 
%    Text for example...
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: 

% Author:       Matthias Althoff, Niklas Kochdumper
% Written:      21-August-2012
% Last update:  29-January-2018
% Last revision:---

%------------- BEGIN CODE --------------

%linearization point p.u of the input is the center of the input u
p.u = ustar;

%linearization point p.x and p.y
p.x = xstar;

%substitute p into the system equation in order to obtain the constant
%input
f0 = sysDisc.mFile(p.x, p.u);

%get jacobian matrices
[A_lin,B_lin] = sysDisc.jacobian(p.x, p.u);

%---------- Made this part as a comment-------------
% uTrans = f0; %B*Ucenter from linOptions.U not added as the system is linearized around center(U)
% %Udelta = B_lin*(U+(-center(U)));
% Udelta = B_lin*(U+(-ustar));
% Uptrans = Udelta + uTrans;
%---------- Made this part as a comment-------------


%------------- END OF CODE --------------