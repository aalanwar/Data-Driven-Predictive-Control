
% P = polytope(Z);
% H = get(P,'H');
% K = get(P,'K');
rand('seed',4500);

clear all
close all
%% system dynamics
%     P = polytope(Z);
%     H = get(P,'H');
%     K = get(P,'K');
dim_x = 5;
A = [-1 -4 0 0 0; 4 -1 0 0 0; 0 0 -3 1 0; 0 0 -1 -3 0; 0 0 0 0 -2];
B_ss = ones(5,1);
C = [1,0,0,0,0];
D = 0;
% define continuous time system
sys_c = ss(A,B_ss,C,D);
% convert to discrete system
samplingtime = 0.05;
sys_d = c2d(sys_c,samplingtime);
initpoints =1;
steps = 300;
totalsamples = initpoints*steps;
%% initial set and input
uref = 8;    
ref = inv(eye(5)-sys_d.A)*sys_d.B*uref;


%X0 = zonotope([ones(dim_x,1),0.1*diag(ones(dim_x,1))]);
%U = zonotope([8,0.0001]);
X0 = zonotope([ref-2,25*diag(ones(dim_x,1))]);
y0 = randPoint(X0);
U = zonotope([uref-1,20-1]);

%noise zontope W
%W = zonotope([zeros(dim_x,1),0.0001*ones(dim_x,1)]);
wfac=0.01;
W = zonotope([zeros(dim_x,1),wfac*ones(dim_x,1)]);

for i=1:size(W.generators,2)
    vec=W.Z(:,i+1);
    GW{i}= [ vec,zeros(dim_x,totalsamples-1)];
    for j=1:totalsamples-1
        GW{j+i}= [GW{i+j-1}(:,2:end) GW{i+j-1}(:,1)];
    end
end

vfac = 0.002;
V = zonotope([zeros(dim_x,1),vfac*ones(dim_x,1)]);
CV = zeros(dim_x,totalsamples);
for i=1:size(V.generators,2)
    vec=V.Z(:,i+1);
     GV{i}= [ vec,zeros(dim_x,totalsamples-1)];
    for j=1:totalsamples-1
        GV{j+i}= [GV{i+j-1}(:,2:end) GV{i+j-1}(:,1)];
    end
end

Vmatzono= matZonotope(CV,GV);
AV = sys_d.A*V;
VAmatzono = sys_d.A*Vmatzono;

% concatinate W to get matrix zonotope Wmatzono
% for i=1:size(W.generators,2)
%     GW1{i}= repmat(W.Z(:,i+1),1,totalsamples);
% end

% for i=1:size(W.generators,2)
%     vec=W.Z(:,i+1);
%     for j=0:totalsamples-1
%         GW{j+i}= [ zeros(dim_x,j),vec,zeros(dim_x,totalsamples-j-1)];
%     end
% end

Wmatzono= matZonotope(zeros(dim_x,totalsamples),GW);


% randomly choose constant inputs for each step / sampling time
for i=1:totalsamples
    u(i) = randPoint(U);
end


% %% simulate the discrete system starting from x0
% x0 = ones(5,1);
% x(:,1) = x0;
% for i=1:totalsamples
%     x(:,i+1) = sys_d.A*x(:,i) + sys_d.B*u(i) +randPoint(W);
% end



x0 = X0.center;
x(:,1) = x0;
index=1;
for j=1:dim_x:initpoints*dim_x
    x(j:j+dim_x-1,1) = randPoint(X0);
    x_v(j:j+dim_x-1,1) =  x(j:j+dim_x-1,1) + randPoint(V);

    for i=1:steps
        utraj(j,i) = u(index);
        x(j:j+dim_x-1,i+1) = sys_d.A*x(j:j+dim_x-1,i) + sys_d.B*u(index) + randPoint(W);
        x_v(j:j+dim_x-1,i+1) =  x(j:j+dim_x-1,i+1) + randPoint(V);
        index=index+1;
    end
end



index_0 =1;
index_1 =1;
for j=1:dim_x:initpoints*dim_x
    for i=2:steps+1
         x_meas_vec_1_v(:,index_1) = x_v(j:j+dim_x-1,i);
        x_meas_vec_1(:,index_1) = x(j:j+dim_x-1,i);
        index_1 = index_1 +1;
    end
    for i=1:steps
        u_mean_vec_0(:,index_0) = utraj(j,i);
        x_meas_vec_0(:,index_0) = x(j:j+dim_x-1,i);
        x_meas_vec_0_v(:,index_0) = x_v(j:j+dim_x-1,i);
        index_0 = index_0 +1;
    end
end

U_full = u_mean_vec_0(:,1:totalsamples); %same as u
Y_0T = x_meas_vec_0_v(:,1:totalsamples);
Y_1T = x_meas_vec_1_v(:,1:totalsamples);
NN = 30;%prediction horizon
Tini = 1;%initial condition horizon
%Hankel matrices
m=1;
p=5;
H_u = data2hankel(U_full,Tini+NN);
H_y = data2hankel(Y_0T,Tini+NN);
U_p = H_u(1:m*Tini,:);
U_f = H_u(m*Tini+1:end,:);
Y_p = H_y(1:p*Tini,:);
Y_f = H_y(p*Tini+1:end,:);

lambda_g = 30;
lambda = 100000;


%% Compute MPC problem
%Cost matrices
N = 6;
Qy = 1e3*eye(5); %define output cost matrix
Qu = 0.001*eye(1);%control cost matrix


t = 1;
% intc = interval([-1;2;1;-1;3],[2;4;4;3;6]);
% consSet = zonotope(intc);
%consSet = zonotope(interval([-10;-10;-10;-10;-10],[10;10;10;10;10]));

%y0 = [-1.9;2.55;3.5;1.9;4.3];

maxsteps = 100;
chosedtimestep = 10;
for timesteps = 1:maxsteps
    if timesteps == 1
        y_t(:,timesteps) = y0;
        YPred(:,1) = y0;
    end
%     else
% 
%     end
    
    
    
    % Control
    uini = sdpvar(m*Tini, 1);
    yini = sdpvar(p*Tini, 1);
    sigma = sdpvar(p*Tini,1);
    g = sdpvar(size(U_p,2),1);
    U_f_times_g = U_f*g;
    Y_f_times_g = Y_f*g;
    u = reshape(U_f_times_g,m,NN); %define control variable based on U_f*g
    y = reshape(Y_f_times_g,p,NN); %define output variable based on Y_f*g

    Constraints = [y_t(:,timesteps) == y(1)];%,...

    
    %for i = 1:N

        Constraints = [Constraints,...
            ([U_p; Y_p]*g == [uini; yini] + [zeros(1*Tini,1);sigma]),...
            % [U_f;Y_f]*g ==[u; y]...
            ];
   % end
    
    %pinv(sys_d.B)*(eye(5)-sys_d.A)*[1;0;0;0;0]

    
    Cost=lambda*norm(sigma,1)+ lambda_g*norm(g,2);
    for i=1:NN
        Cost = Cost + (y(i+1)-ref)'*Qy*(y(i+1)-ref)+ (u(i)-uref)'*Qu*(u(i)-uref);
    end
%    Cost = Cost + (y{i+1}-ref)'*Qy*(y{i+1}-ref);
    options = sdpsettings('verbose',0,'solver','mosek');
    Problem = optimize(Constraints,Cost,options)
    Objective = double(Cost);
    uPred(timesteps) = double(u{1})
    %double(alpha_u(1));
    YPred(:,timesteps+1) = double(y{2});
    %y11(:,timesteps) = double(y{1});
    norm(YPred(:,timesteps)-ref,2);
    xxx=11;
    
    
    Rplotall{timesteps}= zonotope([ double(R{2}.center), double(R{2}.generators)]);
    %%  ploting
    %     projectedDims = {[1 2],[3 4],[4 5]};
    if chosedtimestep == timesteps
        for i =1:N+1
            RoverN{i}= zonotope([ double(R{i}.center), double(R{i}.generators)]) ;
            RoverN_int{i} = interval(RoverN{i});
            yoverN{i} =double(y{i});
            if i<N+1             
                uoverN{i} =double(u{i});
            end
        end
    end
%     for plotRun=1:length(projectedDims)
%         figure; hold on;
%         for i=1:N+1
%           handR=  plot(Rint{i},projectedDims{plotRun},'r');
%         end
%          handC=plot(consSet,projectedDims{plotRun},'k');
%         legend([handR,handC],'R','constrain')
%         % label plot
%         xlabel(['x_{',num2str(projectedDims{plotRun}(1)),'}']);
%         ylabel(['x_{',num2str(projectedDims{plotRun}(2)),'}']);
%     end

    uPred_ = double(u{1}) + double(alpha_u(1))*U.generators;
    
    uPred_final(timesteps) = uPred_;  
    w_point = randPoint(W);
    v_point = randPoint(V);
    y_t(:,timesteps+1) = sys_d.A * y_t(:,timesteps) + sys_d.B * uPred(timesteps) + w_point +v_point - sys_d.A *v_point;
    yt2ref(timesteps)= norm(y_t(:,timesteps)-ref,2)
    halt = 1;
end

projectedDims = {[1 2],[3 4],[4 5]};
for plotRun=1:length(projectedDims)
    figure('Renderer', 'painters', 'Position', [10 10 700 900]); 
    hold on;
    box on;
    for i=1:maxsteps
        handR=  plot(Rplotall{i},projectedDims{plotRun},'r');
        % plot(YPred(projectedDims{plotRun}(1),i),YPred(projectedDims{plotRun}(2),i),'*-k');
    end    
   for i=2:maxsteps+1
            haldy_t= plot(y_t(projectedDims{plotRun}(1),i),y_t(projectedDims{plotRun}(2),i),'+b');
%         else
             haldy_Pred= plot(YPred(projectedDims{plotRun}(1),i),YPred(projectedDims{plotRun}(2),i),'*k');
%         end
    end
    % handC=plot(consSet,projectedDims{plotRun},'k');
    legend([handR,haldy_t,haldy_Pred],'Reachable set $\mathcal{R}_k$','System trajectory $y(k)$','$y$-pred$(k)$','Interpreter','latex')
    % label plot
    xlabel(['x_{',num2str(projectedDims{plotRun}(1)),'}']);
    ylabel(['x_{',num2str(projectedDims{plotRun}(2)),'}']);
    
    
    ax = gca;
    ax.FontSize = 16;
    %set(gcf, 'Position',  [50, 50, 800, 400])
    ax = gca;
    outerpos = ax.OuterPosition;
    ti = ax.TightInset;
    left = outerpos(1) + ti(1);
    bottom = outerpos(2) + ti(2);
    ax_width = outerpos(3) - ti(1) - ti(3)-0.04;
    ax_height = outerpos(4) - ti(2) - ti(4);
    ax.Position = [left bottom ax_width ax_height];
end


%sys_d.A*ref+sys_d.B*double(u{1})-ref
%pinv(sys_d.B)*(eye(5)-sys_d.A)*[1;0;0;0;0]

figure
box on;
han_pred=plot(uPred,'b*-');
hold on
han_uref=plot(uref*ones(size(uPred)),'k-')
ax = gca;
xlabel('Time step $k$','Interpreter','latex')
legend([han_pred,han_uref],'Control input $u^*(k)$','Reference input $r_u(k)$','Interpreter','latex')
ax.FontSize = 16;
%set(gcf, 'Position',  [50, 50, 800, 400])
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3)-0.01;
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
%%-----------------------R over N-------------------------------------------%%

for plotRun=1:length(projectedDims)
    figure('Renderer', 'painters', 'Position', [10 10 700 900]); 
    hold on;
    box on;
    for i=2:N+1
        handR=  plot(RoverN_int{i},projectedDims{plotRun},'r');
        % plot(YPred(projectedDims{plotRun}(1),i),YPred(projectedDims{plotRun}(2),i),'*-k');
    end    
   for i=2:N+1
            haldy_t= plot(yoverN{i}(projectedDims{plotRun}(1)),yoverN{i}(projectedDims{plotRun}(2)),'*k');
%         else
%             haldy_t= plot(y0(projectedDims{plotRun}(1)),y0(projectedDims{plotRun}(2)),'*-k');
%         end
    end
    % handC=plot(consSet,projectedDims{plotRun},'k');
    legend([handR,haldy_t],'Reachable set $\mathcal{R}_k$','$y$-pred$(k)$','Interpreter','latex')
    % label plot
    xlabel(['x_{',num2str(projectedDims{plotRun}(1)),'}']);
    ylabel(['x_{',num2str(projectedDims{plotRun}(2)),'}']);
    
    
    ax = gca;
    ax.FontSize = 16;
    %set(gcf, 'Position',  [50, 50, 800, 400])
    ax = gca;
    outerpos = ax.OuterPosition;
    ti = ax.TightInset;
    left = outerpos(1) + ti(1);
    bottom = outerpos(2) + ti(2);
    ax_width = outerpos(3) - ti(1) - ti(3)-0.01;
    ax_height = outerpos(4) - ti(2) - ti(4);
    ax.Position = [left bottom ax_width ax_height];
end
%%---------------------- plot y's ----------------------------------------%%
figure 
hold on
box on;
%for i =1:dim_x
 handy{1}= plot(y_t(1,:)','r');
 handy_pred{1}= plot(YPred(1,:)','-*r');
 handy{2}= plot(y_t(2,:)','k');
 handy_pred{2}= plot(YPred(2,:)','-*k');
 handy{3}= plot(y_t(3,:)','b');
 handy_pred{3}= plot(YPred(3,:)','-*b');
 handy{4}= plot(y_t(4,:)','g');
 handy_pred{4}= plot(YPred(4,:)','-*g');
 handy{5}= plot(y_t(5,:)','m');
 handy_pred{5}= plot(YPred(5,:)','m-*');
%end

% for i =1:dim_x
%  handy_pred{i}= plot(YPred(i,:)','-*');
% end
xlabel('Time step $k$','Interpreter','latex')
legend([handy{1},handy_pred{1},handy{2},handy_pred{2},handy{3},handy_pred{3},handy{4},handy_pred{4},handy{5},handy_pred{5}],...
        '$y_1(k)$','$y_1$-pred$(k)$','$y_2(k)$','$y_2$-pred$(k)$','$y_3(k)$','$y_3$-pred$(k)$','$y_4(k)$','$y_4$-pred$(k)$','$y_5(k)$','$y_5$-pred$(k)$','Interpreter','latex');
    ax = gca;
    ax.FontSize = 12;
    %set(gcf, 'Position',  [50, 50, 800, 400])
    ax = gca;
    outerpos = ax.OuterPosition;
    ti = ax.TightInset;
    left = outerpos(1) + ti(1);
    bottom = outerpos(2) + ti(2);
    ax_width = outerpos(3) - ti(1) - ti(3)-0.01;
    ax_height = outerpos(4) - ti(2) - ti(4);
    ax.Position = [left bottom ax_width ax_height];
    
    %legend([handy{1},handy{2},handy{3},handy{4},handy{5}],...
%        'y1','y2','y3','y4','y5','Location','northwest');
    
%%-------------------------------------------------------------------- 
figure 
hold on
box on;
plot(yt2ref,'b*-')
xlabel('Time step $k$','Interpreter','latex')
legend('$|| y(k) - r_y(k) ||$','Interpreter','latex')
    ax = gca;
    ax.FontSize = 12;
    %set(gcf, 'Position',  [50, 50, 800, 400])
    ax = gca;
    outerpos = ax.OuterPosition;
    ti = ax.TightInset;
    left = outerpos(1) + ti(1);
    bottom = outerpos(2) + ti(2);
    ax_width = outerpos(3) - ti(1) - ti(3)-0.01;
    ax_height = outerpos(4) - ti(2) - ti(4);
    ax.Position = [left bottom ax_width ax_height];
%%---------------------------------------------------------------------
% % set number of steps in analysis
% totalsteps = 8;
% Y_model = cell(totalsteps+1,1);
% Y_data = cell(totalsteps+1,1);
% % init sets for loop
% Y_model{1} = X0; Y_data{1} = X0;
% 
% for i=1:totalsteps
%     
%     % 1) model-based computation
%     Y_model{i,1}=reduce(Y_model{i,1},'girard',10);
%     Y_model{i+1,1} = sys_d.A * Y_model{i} + sys_d.B * U+W +V + -1*sys_d.A * V;
%     % 2) Data Driven approach
%     Y_data{i,1}=reduce(Y_data{i,1},'girard',400);
%     Y_data{i+1,1} = AB * (cartesianProduct(Y_data{i},U)) +W +V+ -1*AV;
%     
%     
% end
% 
% 
% 
% 
% %% visualization
% 
% projectedDims = {[1 2],[3 4],[4 5]};
% axx{1} = [0.75,1.5,0.5,4]; axx{2} = [0.75,3,0.8,2.2];axx{3} = [0.75,2.3,0.75,2.8];
% index=1;
% numberofplots = 5;%length(X_model)
% for plotRun=1:length(projectedDims)
%     
%     figure('Renderer', 'painters', 'Position', [10 10 700 900])
%     
%     % set axis
%     
%     index=index+1;
%     % plot initial set
%     handleX0 = plot(X0,projectedDims{plotRun},'k-','LineWidth',2);
%     hold on;
%     
%     
%     % plot reachable sets starting from index 2, since index 1 = X0
%     
%     % plot reachable sets from model
%     for iSet=2:numberofplots
%         handleModel=  plot(Y_model{iSet},projectedDims{plotRun},'r');
%     end
%     
%     % plot reachable sets from data
%     for iSet=2:numberofplots
%         handleData=   plot(Y_data{iSet},projectedDims{plotRun},'k');
%     end
%     
%     % label plot
%     xlabel(['x_{',num2str(projectedDims{plotRun}(1)),'}']);
%     ylabel(['x_{',num2str(projectedDims{plotRun}(2)),'}']);
%     %axis(axx{plotRun});
%     % skip warning for extra legend entries
%     warOrig = warning; warning('off','all');
%     legend([handleX0,handleModel,handleData],...
%         'Initial Set','Set from Model','Set from Data','Location','northwest');
%     warning(warOrig);
%     ax = gca;
%     ax.FontSize = 22;
%     %set(gcf, 'Position',  [50, 50, 800, 400])
%     ax = gca;
%     outerpos = ax.OuterPosition;
%     ti = ax.TightInset;
%     left = outerpos(1) + ti(1);
%     bottom = outerpos(2) + ti(2);
%     ax_width = outerpos(3) - ti(1) - ti(3)-0.01;
%     ax_height = outerpos(4) - ti(2) - ti(4);
%     ax.Position = [left bottom ax_width ax_height];
%     %set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
% end
% 



%------------- END OF CODE --------------