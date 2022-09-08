% ZPC: run RMPC with model using polytopes
% 
%
% Inputs:
%    none
%
% Outputs:
%    saved workspace
%
% Example: 
%
% See also: ---

% Author:       Yvonne Stürz 
% Written:      25-March-2021 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------


clear all
close all
rand('seed',4500);

% The system description is xk+1=Axk+Buk+Ewk,yk=Cxk
A = [-1 -4 0 0 0; 4 -1 0 0 0; 0 0 -3 1 0; 0 0 -1 -3 0; 0 0 0 0 -2];
B_ss = ones(5,1);
C = [1,0,0,0,0];
D = 0;
% define continuous time system
sys_c = ss(A,B_ss,C,D);
% convert to discrete system
samplingtime = 0.05;
sys_d = c2d(sys_c,samplingtime);
A = sys_d.A;
B = sys_d.B;

dim_x = 5;
C = eye(5);
E = eye(5);

wfac=0.01;
vfac = 0.002;
Qy = 1e3*eye(5); %define output cost matrix
Qu = 0.001*eye(1);%control cost matrix

uref = 8;    
ref = inv(eye(5)-sys_d.A)*sys_d.B*uref;

y_lb = [-10;2;-10;-10;-10]; %[-2;2;1;-2.5;3];
y_ub = [10;10;10;10;10]; %[2;4;4;3;6];
y0 = [-2;4;3;-2.5;5.5];


% Open loop mini-max solution

N = 2;
U = sdpvar(N,1);
W = sdpvar(N,5);
V = sdpvar(N,5);
x = sdpvar(5,1);

Vw = wfac*ones(5,1);
Vw = [Vw'; -Vw'];
Pw = Polyhedron('V', Vw);
%plot(P1)
Pw.minHRep;
Pw.H;
Pw.He;

test = [0;0;0;0;0];
W_Constraints = [Pw.A * test <= Pw.b; Pw.Ae * test == Pw.be]


Vv = vfac*ones(5,1);
Vv = [Vv'; -Vv'];
Pv = Polyhedron('V', Vv);
%plot(P1)
Pv.minHRep;
Pv.H;
Pv.He;

test = [0;0;0;0;0];
V_Constraints = [Pv.A * test <= Pv.b; Pv.Ae * test == Pv.be]


Y = [];
xk = x;
for k = 1:N
 xk = A*xk + B*U(k)+E*W(k,:)';
 Y = [Y; C*xk + V(k,:)'];
end

F = [kron(ones(N,1),y_lb) <= Y <= kron(ones(N,1),y_ub), kron(ones(N,1),-32) <= U <= kron(ones(N,1),46)];
objective = norm(Y-kron(ones(N,1),ref),2)*Qy(1) + norm(U-kron(ones(N,1),uref),2)*Qu(1);

G = []; 
for k = 1:N
    G = [G, blkdiag(Pw.A,Pv.A) * [W(k,:)'; V(k,:)'] <= [Pw.b; Pv.b], blkdiag(Pw.Ae,Pv.Ae) * [W(k,:)'; V(k,:)'] == [Pw.be;Pv.be]];
end


[Frobust,h] = robustify(F + G,objective,[],[W;V]);


xk = y0;
uk = [];
Y = y0;
ops = sdpsettings;
maxsteps = 80;
execTime=[];
for i = 1:maxsteps
    tic
    optimize([Frobust, x == xk(:,end)],h,ops);
    execTime=[execTime toc];
    xk = [xk A*xk(:,end) + B*value(U(1)) + E*wfac*(-1+2*rand(1)*ones(5,1))];
    Y = [Y, C*xk(:,end) + vfac*(-1+2*rand(1)*ones(5,1))];
    uk = [uk value(U(1))];
end

Cost_rob_ol_tot=0;
Cost_rob_ol=[];
for i = 1:maxsteps
    Cost_rob_ol = [Cost_rob_ol, (Y(:,i+1)-ref)'*Qy*(Y(:,i+1)-ref)+ (uk(:,i)-uref)'*Qu*(uk(:,i)-uref)];
    Cost_rob_ol_tot = Cost_rob_ol_tot + (Y(:,i+1)-ref)'*Qy*(Y(:,i+1)-ref)+ (uk(:,i)-uref)'*Qu*(uk(:,i)-uref);
    yt2ref_poly(i) = norm(Y(:,i)-ref,2);
end
Cost_rob_ol_tot
meanExecTime=mean(execTime)
stdExecTime= std(execTime)
% figure(1)
% %plot([C*xk  + vfac*(-1+2*rand(1)*ones(5,1))]')
% plot([Y]')
% hold on, plot(kron(ones(1,100),ref)')
% figure(2)
% %plot([C*xk  + vfac*(-1+2*rand(1)*ones(5,1))]')
% plot([Y]')
% hold on, plot(kron(ones(1,100),ref)')
% hold on, plot(kron(ones(1,100),y_lb)') 
% hold on, plot(kron(ones(1,100),y_ub)') 

save('workspaces\poly');

