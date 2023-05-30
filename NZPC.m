% NZPC: run NZPC with model and without the model
% i.e., run: 
% 1- Nonlinear robust data driven Predictive Control scheme (NZPC) 
% 2-Same NZPC while knowing the model (RMPC-zono)

% Author:       Mahsa Farjadnia, Amr Alanwar
% Last update:  20-May-2023

%------------- BEGIN CODE --------------

rand('seed',45);
%% 
clear all;
close all;
clc;
%% 
%----MPC horizon
N = 3; 
%----Simulation steps 
maxsteps = 150;
%% Flag to include terminal set 
terminalSetFlag = 0;
%% Initial Zonotopes to generate data
%----Output system matrix
H  = [1 0.001;-0.01 1];
%----Initial conditions
x0 = [-2;-20.5];
y0 = H * x0;
%----Initial state zonotope
X0 = zonotope([x0,diag([0.01;1])]);

%----Input zonotope

uref  = [0;0.007];
u_gnr = diag([5;3]); 
U     = zonotope([uref,u_gnr]); 

%----Dimensions of x,y,u
dim_x = size(x0,1);
dim_y = size(y0,1);
dim_u = size(uref,1);


%----Define output cost matrix
Qy = 5*eye(dim_y); 
%----Control cost matrix
Qu = 0.02*eye(dim_u);
%----Output Constraints

y_lb = [-3;-22];  
y_ub = [0.25;2.7];
intc = interval(y_lb,y_ub);

%----Reference output
yref = [0;0];

%----Number of initial states
initpoints = 50;
%----Number of steps for each trajectory
steps = 10;
%----Total number of samples
totalsamples = initpoints*steps;
%% Define System
%----CSTR system in cont time
dt  = 0.015;                    %sample time
fun = @(x,u) cstrDiscr(x,u,dt); %discrete time function

options_MD.tensorOrder = 2;
options_MD.errorOrder  = 5;

%----Compute discretized CSTR system 
sysDisc = nonlinearDT('stirredTankReactor',fun,dt,dim_x,dim_u);

[A_lin,B_lin] = sysDisc.jacobian(x0, uref); 

%% Other zonotopes
%----Inf norm of the states
eta   = 22;
R_eta =  zonotope([[0;0],eta*eye(dim_x)]);
%----Zonotope of the Process noise (w)

% Generator related factor (zonotope)
wFac = [2e-4;2e-2]; 

% Center of noise factor (zonotope)
wCenFac = 0; 

% Noise zonotope 
W  = zonotope([wCenFac*ones(dim_x,1),wFac]);

% Center of matrix zonotope 
CW = wCenFac*ones(dim_x,totalsamples);

for i = 1 : size(W.generators,2)
    vec   = W.Z(:,i+1);
    GW{i} = [ vec,zeros(dim_x,totalsamples-1)];
    for j = 1 : totalsamples-1
        GW{j+i} = [GW{i+j-1}(:,2:end) GW{i+j-1}(:,1)];
    end
end

%----Matrix zonotpe of noise w (M_w)
Wmatzono = matZonotope(CW,GW);

%----Zonotope of measurement noise (V)

% Generator related factor (zonotope)
vFac    = [1e-3;1e-2];
% Center of noise factor (zonotope)
vCenFac = 0;

% Noise zonotope 
V  = zonotope([vCenFac*ones(dim_x,1),vFac]);

% Center of matrix zonotope 
CV = vCenFac*ones(dim_x,totalsamples);

for i = 1 : size(V.generators,2)
    vec   = V.Z(:,i+1);
    GV{i} = [ vec,zeros(dim_x,totalsamples-1)];
    for j = 1:totalsamples-1
        GV{j+i} = [GV{i+j-1}(:,2:end) GV{i+j-1}(:,1)];
    end
end
%----Matrix zonotpe of measurement noise v (M_v)
Vmatzono = matZonotope(CV,GV);

%% Generate Data
%----Generate u
for i = 1 : totalsamples
    u(:,i) = randPoint(U);
end

%----Get state and output trajectories
x(:,1) = randPoint(X0);
y(:,1) = H*x(:,1);

index = 1;
for j = 1 : dim_x:initpoints*dim_x
    x(j:j+dim_x-1,1) = randPoint(X0);
    y(j:j+dim_x-1,1) = H*x(j:j+dim_x-1,1) + randPoint(V);
    for i = 1 : steps
        x(j:j+dim_x-1,i+1) = fun(x(j:j+dim_x-1,i),u(:,index)) + randPoint(W);
        y(j:j+dim_x-1,i+1) =  H*x(j:j+dim_x-1,i+1) + randPoint(V);
        index = index+1;
    end
end


%----Combine trajectories (Y+,Y-,X+,X-)
index_0 = 1;
index_1 = 1;

for j = 1 : dim_x:initpoints*dim_x
    for p = 1 : dim_y:initpoints*dim_y
        for i = 2 : steps+1
            x_meas_vec_1(:,index_1) = x(j:j+dim_x-1,i);
            y_meas_vec_1(:,index_1) = y(p:p+dim_y-1,i);
            index_1 = index_1 +1;
        end
        for i = 1 : steps
            x_meas_vec_0(:,index_0) = x(j:j+dim_x-1,i);
            y_meas_vec_0(:,index_0) = y(p:p+dim_y-1,i);
            index_0 = index_0 +1;
        end
    end
end
% U_data is U_-, Y_0T is Y_- , Y_1T is Y_+
U_data = u(:,1:totalsamples); %same as u
Y_0T   = y_meas_vec_0(:,1:totalsamples);
Y_1T   = y_meas_vec_1(:,1:totalsamples);
%% Compute Lipschitz constants
ZepsFlag = 1;
stepsLip = 1;
initpointsLip = 200;
[gamma,L] = compLipConst(fun,U,X0,stepsLip,initpointsLip,dim_x);

eps(1) = L(1) .* gamma(1)/2;
eps(2) = L(2) .* gamma(2)/2;

Zeps   = zonotope([zeros(2,1),diag(eps)])

%% Compute NZPC problem

%----Data-driven-----%
x_t(:,1)   = x0;
y_t(:,1)   = y0;
YPred(:,1) = y0;
uPred(:,1) = uref;
%----model-based-----% 
x_t_model(:,1)   = x0;
y_t_model(:,1)   = y0;
YPred_model(:,1) = y0;
uPred_model(:,1) = uref;
%-------Noise--------%
w_point = randPoint(W);
v_point = randPoint(V);

timesteps = 1;
%----data-driven-----%
x_t(:,timesteps+1)   = fun(x_t(:,timesteps),uPred(:,timesteps)) + w_point; 
y_t(:,timesteps+1)   = H*x_t(:,timesteps+1) + v_point; 
Rplotall{timesteps}  = interval(y_t(:,timesteps+1));
YPred(:,timesteps+1) = y_t(:,timesteps+1);
yt2ref(timesteps)    = norm(y_t(:,timesteps)-yref,2);

%----model-based-----%
x_t_model(:,timesteps+1)   = fun(x_t_model(:,timesteps),uPred_model(:,timesteps)) + w_point; 
y_t_model(:,timesteps+1)   = H*x_t_model(:,timesteps+1) + v_point;
Rplotall_model{timesteps}  = interval(y_t_model(:,timesteps+1));
YPred_model(:,timesteps+1) = y_t_model(:,timesteps+1);
yt2ref_model(timesteps)    = norm(y_t_model(:,timesteps)-yref,2);

%----Required parameters for computing Reachable sets-----%

%H_pinv = H'inv(HH') output matrix
H_pinv = H'*inv(H*H');
%g_eta = I_nx -H_pinv*H 
g_eta   = eye(dim_x)-H_pinv*H;

for timesteps = 2 : maxsteps
    
    %% NZPC  NZPC NZPC NZPC NZPC NZPC NZPC NZPC NZPC NZPC NZPC NZPC

    U_data = [U_data(:,2:end) uPred(:,timesteps-1)]; 
    Y_0T   = [Y_0T(:,2:end) y_t(:,timesteps-1)];
    Y_1T   = [Y_1T(:,2:end) y_t(:,timesteps)];    
    %----Linearization point
    yStar  = y_t(:,timesteps-1);   
    uStar  = U.center;
    vStar  = vCenFac*ones(dim_x,1); 
    xqstar = H_pinv*(yStar-vStar) + g_eta*R_eta.center;
    yStarMat  = repmat(yStar,1,size(Y_0T,2));
    uStarMat  = repmat(uStar,1,size(U_data,2));
    vStarMat  = repmat(vStar,1,size(U_data,2));
    xqStarMat = repmat(xqstar,1,size(U_data,2));
    %----Computing approximate linearized model
    oneMat = repmat([1],1,size(U_data,2));
    y_tmp  = H_pinv*Y_0T + (-1)*xqStarMat;
    Mprime = (Y_1T - H *Wmatzono.center )*pinv([oneMat; y_tmp ;U_data+(-1)*uStarMat]);
    %----Computing model mismatch & the Lagrange remainder
    
    min_l = Y_1T(:,1)- Mprime * [0;H_pinv*Y_0T(:,1)-xqstar;U_data(:,1)+-1*uStar];
    max_l = Y_1T(:,1)- Mprime * [0;H_pinv*Y_0T(:,1)-xqstar;U_data(:,1)+-1*uStar];
 
    for i = 1 : size(Y_1T,1)
        for j = 2 : size(Y_1T,2)
            tmp_var_0 = Mprime * [0;H_pinv*Y_0T(:,j)-xqstar;U_data(:,j)+-1*uStar];
            tmp_var(i,j) = Y_1T(i,j)- tmp_var_0(i,:);
             if tmp_var(i,j) < min_l(i,1)
                min_l(i,1) = tmp_var(i,j);
            elseif tmp_var(i,j) > max_l(i,1)
                max_l(i,1) = tmp_var(i,j);
             end
        end
    end
    cent_zon_noise = [1; -H_pinv*V.center + g_eta*R_eta.center; zeros(size(U_data,1),1)];
    gen_zon_noise_1 = [zeros(1,size(V.generators,2)); -H_pinv*V.generators; zeros(size(U_data,1),size(V.generators,2))];
    gen_zon_noise_2 = [zeros(1,size(R_eta.generators,2)); g_eta*R_eta.generators; zeros(size(U_data,1),size(R_eta.generators,2))];
    noise_zonotope = zonotope([cent_zon_noise,[gen_zon_noise_1,gen_zon_noise_2]]);
    Z_L = zonotope(interval(min_l,max_l))+ -Mprime*noise_zonotope+ -H_pinv*W+ -1*V;

    %% Solving the optimization problem
    %----Define sdpvar variables
    u = sdpvar(dim_u*ones(1,N),ones(1,N));
    y = sdpvar(dim_y*ones(1,N+1),ones(1,N+1));
    alpha_u_1 = sdpvar(1,N);
    alpha_u_2 = sdpvar(1,N);

    sinf = sdpvar(dim_y*ones(1,N+1),ones(1,N+1));
    ssup = sdpvar(dim_y*ones(1,N+1),ones(1,N+1));
    R = {};
    R{1} = zonotope([y_t(:,timesteps)]);
    
    %----Define Constaints
    
    % Set the first constraint as y_t = current measured y
    Constraints = [y{1} == y_t(:,timesteps) ];
          
    for i = 1:N    
        U_m_ustar = zonotope(u{i}-uStar);  
        Z_temp    = cartesianProduct( zonotope(H*R{i}.center+(-1)*H*V.center +g_eta *R_eta.center...
                    +(-1)*xqstar,[H*R{i}.generators,H*V.generators, g_eta *R_eta.generators]),U_m_ustar );
        
        Mprime_mul_zono = Mprime*cartesianProduct(zonotope([1]), Z_temp);
        noiseLipZono    = Z_L + H*W + V;
        R{i+1} = zonotope([Mprime_mul_zono.center+noiseLipZono.center,[noiseLipZono.generators,Mprime_mul_zono.generators]]);
    end
    
    for i = 2:N+1 
       R{i} = R{i} + Zeps;
    end
    
    for i = 1:N
        RZ = R{i+1}.Z;
        c  = RZ(:,1);
        
        %Determine left and right limit of the reahable set (convert to interval)
        delta = sum(abs(RZ),2) - abs(c);
        leftLimit_zpc{i}  = c - delta;
        rightLimit_zpc{i} = c + delta;
        
        %Specify all constraints
        Constraints = [Constraints,...
            u{i} == U.center + alpha_u_1(i) * u_gnr(:,1) + alpha_u_2(i) * u_gnr(:,2),...
            leftLimit_zpc{i}  >= intc.inf,...
            rightLimit_zpc{i} <= intc.sup,...
            y{i+1}  >= leftLimit_zpc{i},...
            y{i+1}  <= rightLimit_zpc{i},...
            alpha_u_1(i) <= 1 , ...
            alpha_u_1(i) >= -1, ...
            alpha_u_2(i) <= 1 , ...
            alpha_u_2(i) >= -1, ...
            ];
    end

    %----Difine the objective function of NZPC
    Cost = 0;
    for i = 1 : N
        Cost = Cost + (y{i+1}-yref)'*Qy*(y{i+1}-yref)+ (u{i}-uref)'*Qu*(u{i}-uref);
    end   
    %----Solve NZPC
    Problem   = optimize(Constraints,Cost)
    Objective = double(Cost);
    Status_ZPC{timesteps} = Problem.info;
    uPred(:,timesteps)    = double(u{1});
    YPred(:,timesteps+1)  = double(y{2});

    %----Save reachable sets for plotting
    Rplotall{timesteps} = interval(zonotope([ double(R{2}.center), double(R{2}.generators)]));    
    %% NZPC given the model:RMPC-zono RMPC-zono RMPC-zono RMPC-zono RMP
    Rx = {};
	R  = {};
    Rx{1} = zonotope([x_t_model(:,timesteps)]);
    R{1}  = zonotope([y_t_model(:,timesteps)]);
    
    %----Linearization point
    xStar = x_t_model(:,timesteps-1);
    uStar = U.center;
    vStar = vCenFac*ones(dim_x,1); 
    
    %----Compute the linearized dynamics
    [f0,A_lin,B_lin] = linearize_DTT(sysDisc,xStar,uStar);
    controllability(:,timesteps) = length(A_lin) - rank(ctrb(A_lin,B_lin));
    observability(:,timesteps)   = length(A_lin) - rank(obsv(A_lin,H));

    Verror = linError_mixed_noInt_DTT(sysDisc, Rx{1},U,xStar,uStar);
    %% Solve the optimization problem
    %----Define sdpvar variables
    alpha_model_u_1 = sdpvar(1,N);
    alpha_model_u_2 = sdpvar(1,N);
    u_model   = sdpvar(dim_u*ones(1,N),ones(1,N));
    y_model   = sdpvar(dim_y*ones(1,N+1),ones(1,N+1));
    sinf = sdpvar(dim_y*ones(1,N+1),ones(1,N+1));
    ssup = sdpvar(dim_y*ones(1,N+1),ones(1,N+1));
    
    %----Define Constaints
    
    % Set the first constraint as y_t_model = current measured y
    Constraints = [y_model{1} == y_t_model(:,timesteps)];
    %Specify all constraints
    for i = 1:N    
        %----Over-approximating the model-based reachable sets
        
        %Linearized model
        Uptrans = B_lin*(u_model{i}-uStar) + f0;        
        AR_delta    = A_lin*(Rx{i}-xStar);
        Uptrans_znt = zonotope(Uptrans);
        f_lin = zonotope([AR_delta.center + Uptrans_znt.center ,[AR_delta.generators Uptrans_znt.generators]]);
        
        % Over-approximated state reachable set Rx{i}+error+w(i) 
        Rest_R  = zonotope([Verror.center + W.center,[Verror.generators,W.generators]]);
        Rx{i+1} = zonotope([f_lin.center+Rest_R.center,[f_lin.generators,Rest_R.generators]]);
        
        % Over-approximated output reachable set H*Rx + V;
        R{i+1} = zonotope([H*(Rx{i+1}.center+V.center),[H*Rx{i+1}.generators,H*V.generators]]);
        
        % Extract center
        c = R{i+1}.Z(:,1);
        
        % Determine left and right limit
        delta = sum(abs(R{i+1}.Z),2) - abs(c);
        leftLimit_rmpc{i}  = c - delta;
        rightLimit_rmpc{i} = c + delta;
        
        %Specify all constraints
        Constraints = [Constraints,...
            u_model{i} == U.center + alpha_model_u_1(i) * u_gnr(:,1) + alpha_model_u_2(i) * u_gnr(:,2),...           
            leftLimit_rmpc{i}  >= intc.inf,...
            rightLimit_rmpc{i} <= intc.sup,...
            y_model{i+1}  >= leftLimit_rmpc{i},...
            y_model{i+1}  <= rightLimit_rmpc{i},...            
            alpha_model_u_1(i) <= 1 , ...
            alpha_model_u_1(i) >= -1, ...
            alpha_model_u_2(i) <= 1 , ...
            alpha_model_u_2(i) >= -1, ...
            ];
    end
    
    %----Difine the objective function of RMPC-zono
    Cost_model = 0;
    for i = 1:N
        Cost_model = Cost_model + (y_model{i+1}-yref)'*Qy*(y_model{i+1}-yref) + (u_model{i}-uref)'*Qu*(u_model{i}-uref);
    end  
   %----Solve RMPC-zono optimization problem
    Problem    = optimize(Constraints,Cost_model)
    Objective  = double(Cost_model);
    uPred_model(:,timesteps)   = double(u_model{1});
    YPred_model(:,timesteps+1) = double(y_model{2});
    Status_RMPC{timesteps} = Problem.info;
    %----Save model-based reachable sets for plotting
    Rplotall_model{timesteps}  = interval(zonotope([ double(R{2}.center), double(R{2}.generators)]));   
    
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %----Apply the optimal control input to the plant 
    w_point = randPoint(W);
    v_point = randPoint(V);
    % Consider NZPC approach
    x_t(:,timesteps+1) = fun(x_t(:,timesteps),uPred(:,timesteps)) + w_point; 
    y_t(:,timesteps+1) = H*x_t(:,timesteps+1)+ v_point; 
    % Consider RMPC-zono approach
    x_t_model(:,timesteps+1) = fun(x_t_model(:,timesteps),uPred_model(:,timesteps)) + w_point; 
    y_t_model(:,timesteps+1) = H*x_t_model(:,timesteps+1) + v_point;
    
    % Save norm(y(t)-y_ref)
    yt2ref(timesteps) = norm(y_t(:,timesteps)-yref,2);
    yt2ref_model(timesteps)  = norm(y_t_model(:,timesteps)-yref,2);
end
%----Compute the objective function values for NZPC & RMPC-zono
Cost=0;
for i=1:timesteps
    Cost_vec(i) = (y_t(:,i+1)-yref)'*Qy*(y_t(:,i+1)-yref)+ (uPred(:,i)-uref)'*Qu*(uPred(:,i)-uref);
    Cost = Cost + Cost_vec(i);
end

Cost_model = 0;
for i = 1:timesteps
    
    Cost_model_vec(i) = (y_t_model(:,i+1)-yref)'*Qy*(y_t_model(:,i+1)-yref)+ (uPred_model(:,i)-uref)'*Qu*(uPred_model(:,i)-uref);
    Cost_model = Cost_model + Cost_model_vec(i);
end
%% Additional evaluations
%----Check if all points are within the reachable regions. 
perfectFlagZPC = 1;
for i = 2:maxsteps
   if(~Rplotall{i}.in(y_t(:,i+1)))
       perfectFlagZPC = 0;
       fprintf('Warning: output trajectory is not within the predicted output reachable set.\n');
       break;
   end 
end

if perfectFlagZPC == 1
    fprintf('Output trajectory is within the predicted output reachable set.\n');
end

%----Check for observability and Controllability conditions
observabilityFlag = 1;
for ee = 2:maxsteps
    if observability(:,ee) == 1
        observabilityFlag = 0;
        fprintf('Warning: the linearized model is not observable at time %d.\n',ee);
        break;
    end
end

if observabilityFlag == 1
    fprintf('All linearized models within the simulation time are observable.\n');
end
controllabilityFlag = 1;
for ff = 2:maxsteps
    if controllability(:,ff) == 1
       controllabilityFlag = 0;
       fprintf('Warning: the linearized model is not controllable at time %d.\n',ff);       
        break;
    end
end
if controllabilityFlag == 1
    fprintf('All linearized models within simulation time are controllable.\n');
end
%----Check the feasibility of the optimization problems for mosek
% Status_RMPCFlag = 1;
% for b = 2:maxsteps
%     if ~strcmp(Status_RMPC{b},'Successfully solved (MOSEK)')
%        Status_RMPCFlag = 0;
%        fprintf('Warning: optimization problem (RMPC-zono) did not successfully solved at time %d.\n',b);
%         break;
%     end
% end
% 
% if Status_RMPCFlag == 1
%     fprintf('All optimization problems (RMPC-zono) within the simulation time were successfully solved.\n');
% end
% 
% Status_ZPCFlag = 1;
% for gg = 2:maxsteps
%     if ~strcmp(Status_ZPC{gg},'Successfully solved (MOSEK)')
%        Status_ZPCFlag = 0;
%        fprintf('Warning: optimization problem (NZPC) did not successfully solved at time %d.\n',gg);
%         break;
%     end
% end
% if Status_ZPCFlag == 1
%     fprintf('All optimization problems (NZPC) within the simulation time were successfully solved.\n');
% end
 
%% ----------------- Plot reachable sets over steps--------------------%%
projectedDims = {[1 2]};
for plotRun=1:length(projectedDims)
    figure('Renderer', 'painters', 'Position', [10 10 700 800]);
    hold on;
    box on;
    for i = 2:maxsteps
        handR =  plot(Rplotall{i},projectedDims{plotRun},'r');
    end
    for i = 3:maxsteps+1
        haldy_t = plot(y_t(projectedDims{plotRun}(1),i),y_t(projectedDims{plotRun}(2),i),'+b');
       
        haldy_Pred = plot(YPred(projectedDims{plotRun}(1),i),YPred(projectedDims{plotRun}(2),i),'*k');
    
    end
    warOrig = warning; warning('off','all');
    legend([handR,haldy_t,haldy_Pred],'Reachable set $\hat{\mathcal{R}}_{t+k+1|t}^y$','System trajectory $y(k)$','$y_{pred}(k)$','Interpreter','latex','Location','southeast')
    warning(warOrig);
    xlabel(['$y_{',num2str(projectedDims{plotRun}(1)),'}$'],'Interpreter','latex');
    ylabel(['$y_{',num2str(projectedDims{plotRun}(2)),'}$'],'Interpreter','latex');
    
    
    ax = gca;
    outerpos = ax.OuterPosition;
    ax.FontSize = 19;
    ti = ax.TightInset;
    left = outerpos(1) + ti(1)-0.015;
    bottom = outerpos(2) + ti(2);
    ax_width = outerpos(3) - ti(1) - ti(3)-0.01;
    ax_height = outerpos(4) - ti(2) - ti(4);
    ax.Position = [left bottom ax_width ax_height];
end
ax.YLim = [-18 2.78];
ax.XLim = [-2 0.5];
axes('position',[0.15 0.55 0.3 0.4])
box on 
for plotRun=1:length(projectedDims)
    hold on;
    box on;
    for i = 2:maxsteps
        handR =  plot(Rplotall{i},projectedDims{plotRun},'r');
    end
    for i = 3:maxsteps+1
        haldy_t = plot(y_t(projectedDims{plotRun}(1),i),y_t(projectedDims{plotRun}(2),i),'+b');
       
        haldy_Pred = plot(YPred(projectedDims{plotRun}(1),i),YPred(projectedDims{plotRun}(2),i),'*k');
    
    end
end
ax = gca;

ax.YLim = [-2.2 2.78];
ax.XLim = [-0.13 0.24];
hold on 
annotation('arrow',[0.735238095238095 0.536190476190476],...
    [0.870666666666667 0.831666666666667]);
%% ----------------------------Plot u---------------------------------%%
figure('Renderer', 'painters', 'Position', [10 10 700 800]);
box on;
hold on;
han_pred = plot(uPred(1,:),'b*-');
han_pred_model = plot(uPred_model(1,:),'r+-');
han_uref = plot(uref(1).*ones(size(uPred(1,:))),'k-');
ax = gca;
xlabel('Time step $k$','Interpreter','latex')

legend([han_pred,han_pred_model,han_uref],'NZPC $u_1^*(k)$','RMPC-zono $u_1^*(k)$','Reference input $r_{u_1}(k)$','Interpreter','latex','Location','southeast')
ax.FontSize = 19;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3)-0.01;
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

%-----------------------
figure('Renderer', 'painters', 'Position', [10 10 700 800]);
box on;
hold on
han_pred = plot(uPred(2,:),'b*-');
han_pred_model = plot(uPred_model(2,:),'r+-');
han_uref = plot(uref(2).*ones(size(uPred(2,:))),'k-');
ax = gca;
xlabel('Time step $k$','Interpreter','latex')
legend([han_pred,han_pred_model,han_uref],'NZPC $u_2^*(k)$','RMPC-zono $u_2^*(k)$','Reference input $r_{u_2}(k)$','Interpreter','latex','Location','southeast')
ax.FontSize = 19;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1)+0.01;
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3)-0.01;
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
%% ------------------------- Cost -------------------------------------%%
figure('Renderer', 'painters', 'Position', [10 10 700 800]);
hold on
box on;
han_yt2ref_model=plot(Cost_model_vec,'r+-');
han_yt2ref=plot(Cost_vec,'b*-');
xlabel('Time step $k$','Interpreter','latex')
legend([han_yt2ref,han_yt2ref_model],'NZPC cost','RMPC-zono cost','Interpreter','latex')
ax = gca;
ax.FontSize = 19;
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1)+0.01;
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3)-0.01;
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

%% -------------------------- Plot y_1-------------------------------%%
figure('Renderer', 'painters', 'Position', [10 10 700 800]);
hold on
box on;
hand_y_t= plot(y_t(1,:),'b*-');
hand_y_t_model= plot(y_t_model(1,:),'r+-');
%hand_YPred_t=plot(YPred(1,:),'b+-');
hand_Yref = plot(yref(1).*ones(size(YPred(1,:))),'k-');

handcon = plot(intc.inf(1)*ones(size(y_t(1,:))),'k--');
handcon = plot(intc.sup(1)*ones(size(y_t(1,:))),'k--');
axis([1,maxsteps ,min(min(y_t(1,:)),intc.sup(1))-1, intc.sup(1)+1]);
xlabel('Time step $k$','Interpreter','latex')
%ylabel('$y_1(k)$','Interpreter','latex')

legend([hand_y_t,hand_y_t_model,handcon,hand_Yref],'NZPC $y_1(k)$','RMPC-zono $y_1(k)$','constraint','reference output $y_{r_1}(k)$','Interpreter','latex','Location','southeast');
ax = gca;
ax.FontSize = 19;
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3)-0.01;
ax_height = outerpos(4) - ti(2) - ti(4)-0.01;
ax.Position = [left bottom ax_width ax_height];
ax.YLim = [-2.03 0.3];
ax.XLim = [1 150];
%% --------------------------- Plot y_2------------------------------%%
figure('Renderer', 'painters', 'Position', [10 10 700 800]);
hold on
box on;
hand_y_t= plot(y_t(2,:),'b*-');
hand_y_t_model= plot(y_t_model(2,:),'r+-');
hand_Yref = plot(yref(2).*ones(size(YPred(2,:))),'k-');

handcon = plot(intc.inf(2)*ones(size(y_t(2,:))),'k--');
handcon = plot(intc.sup(2)*ones(size(y_t(2,:))),'k--');
axis([1,maxsteps ,min(min(y_t(2,:)),intc.sup(2))-1, intc.sup(2)+1]);
xlabel('Time step $k$','Interpreter','latex')
legend([hand_y_t,hand_y_t_model,handcon,hand_Yref],'NZPC $y_2(k)$','RMPC-zono $y_2(k)$','constraint','reference output $y_{r_2}(k)$','Interpreter','latex','Location','southeast');
ax = gca;
ax.FontSize = 19;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3)-0.01;
ax_height = outerpos(4) - ti(2) - ti(4)-0.01;
ax.Position = [left bottom ax_width ax_height];

