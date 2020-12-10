
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
%X0 = zonotope([ref-1,0.2*diag(ones(dim_x,1))]);
X0 = zonotope([ref-2,25*diag(ones(dim_x,1))]);
c1 = X0.Z(:,1);

%determine left and right limit
delta1 = sum(abs(X0.Z),2) - abs(c1);
% leftLimit = c - delta;
% rightLimit = c + delta;
intc = interval([-4;2;1;-3;4],[-1;4;4;3;6]);
% intc = interval([-10;1;1;-5;-5],[10;4;5;5;10]);
% intc = interval([-5;-5;-5;-5;-5],[5;5;5;5;5]);
% intc_delta = interval([-2;2;1;-1;3]+delta1,[2;4;4;3;6]-delta1);
%y0 = randPoint(X0&zonotope(intc_delta));
y0 = [-3;3.5;3;-2.5;5.5];
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




% plot simulated trajectory
figure;
subplot(1,2,1); hold on; box on; plot(x(1,:),x(2,:),'b'); xlabel('x_1'); ylabel('x_2');
subplot(1,2,2); hold on; box on; plot(x(3,:),x(4,:),'b'); xlabel('x_3'); ylabel('x_4');
close;


AB = (Y_1T + -1* Vmatzono + -1*Wmatzono+VAmatzono)*pinv([Y_0T;U_full]);
%X1W_cen =  Y_1T - Wmatzono.center;
%X1W = matZonotope(X1W_cen,Wmatzono.generator);
% XU = [X_0T;U_full];
% XU_inv = XU'*inv(XU*XU');

intAB11 = intervalMatrix(AB);
intAB1 = intAB11.int;
intAB1.sup >= [sys_d.A,sys_d.B]
intAB1.inf <= [sys_d.A,sys_d.B]





%% Compute MPC problem
%Cost matrices
N = 3;
Qy = 1e3*eye(5); %define output cost matrix
Qu = 0.001*eye(1);%control cost matrix


t = 1;

% consSet = zonotope(intc);
%consSet = zonotope(interval([-10;-10;-10;-10;-10],[10;10;10;10;10]));

%y0 = [-1.9;2.55;3.5;1.9;4.3];

maxsteps = 40;
for timesteps = 1:maxsteps
    if timesteps == 1
        y_t(:,timesteps) = y0;
        prev_v_point =zeros(5,1);
    end
%     else
% 
%     end
    
    
    
    % Control
    sinf = sdpvar(1);
    ssup = sdpvar(1);
    
    u = sdpvar(1*ones(1,N),ones(1,N));
    y = sdpvar(5*ones(1,N+1),ones(1,N+1));
    zono_order = 20;
    beta = sdpvar((zono_order*5+1)*ones(1,N),ones(1,N));
    beta_u = sdpvar;
    alpha_u = sdpvar*ones(1,N); 
    alpha_w = sdpvar(1,N+1);
    alpha_v = sdpvar(1,N+1);
    alpha_av = sdpvar(1,N+1);
    alpha_w_2 = sdpvar(1,N+1);
    R={};
    R{1} = zonotope([y_t(:,timesteps)]);
    
    Constraints = [y_t(:,timesteps) == y{1}];%,...
                % y{1}+alpha_w_2(1) * W.generators >= intc.inf,...
                % y{1}+alpha_w_2(1) * W.generators <= intc.sup];
    
    %sdpvar w_cen(5,1);
    %sdpvar w_gen(5,1);

    %W_sdp = zonotope([w_cen,w_gen]);
    %0.005*ones(dim_x,1)
    %AB= reduce(AB,'girard',1,0);%0.4,0);
    AB_gen =[];
    
    for i = 1:N
%         if i>1
%             R{i}= and(R{i},consSet,'averaging');
%         end
        %card_cen = [y{i};u{i}];
        card_cen = [R{i}.center;u{i}];
        genleni = size(R{i}.generators,2);
        card_zono = zonotope([card_cen,[R{i}.generators;zeros(1,genleni)]]);
        %R{i+1} = AB * (cartesianProduct(R{i}, zonotope([u{i},0]) )) +W_sdp;
        ABcard = intervalMatrix(AB)* card_zono;
         R{i+1} = zonotope([ABcard.center,[ABcard.generators,W.generators,V.generators,AV.generators]]);%AB * card_zono + W_sdp;
        %convert R to interval
        %extract center
        %R{i+1}= and(R{i+1},consSet,'averaging');
        c = R{i+1}.Z(:,1);
        
        %determine left and right limit
        delta = sum(abs(R{i+1}.Z),2) - abs(c);
        leftLimit{i} = c - delta;
        rightLimit{i} = c + delta;
        
        
        genlenip1 = length(R{i+1}.generators);
        Constraints = [Constraints,...
           %  y{i+1} <= 100,...
          %   y{i+1} >= -100,...
            u{i} == U.center + alpha_u(i) * U.generators,...
       %      u{i} >= 8,...
            % W_sdp.center == zeros(dim_x,1);
            % W_sdp.generators == 0.005*ones(dim_x,1);
            %u{i} == U.center + U.generators*beta_u;
            %beta_u >= -1;
            %beta_u <= 1;
            %y{i+1} == sys_d.A *y{i} + sys_d.B*u{i},...
            %y{i+1} == R{i+1}.center + R{i+1}.generators*beta{i}(1:genlenip1),...
   %         y{i+1} >= leftLimit{i},...
   %         y{i+1} <= rightLimit{i},...
          %  y{i+1}   - 1.1*delta >= intc.inf,...
          %  y{i+1}   + 1.1*delta <= intc.sup,... 
            y{i+1}   - sinf >= intc.inf,...
            y{i+1}   + ssup <= intc.sup,... 
            sinf >= 0,...
            ssup <= 0,...
          %   y{i+1}   - delta >= intc.inf,...
          %  y{i+1}   + delta <= intc.sup,...
           %  y{i+1}   >= intc.inf,...
           % y{i+1}    <= intc.sup,...
           % y{i+1}+ (alpha_w_2(i) * W.generators) + alpha_v(i) * V.generators + alpha_av(i) * AV.generators  >= intc.inf,...
          %  y{i+1}+ (alpha_w_2(i) * W.generators) + alpha_v(i) * V.generators + alpha_av(i) * AV.generators <= intc.sup,...            
           % beta{i}(1:genlenip1) >= -1*ones(genlenip1,1),...
           % beta{i}(1:genlenip1) <= ones(genlenip1,1),...
            alpha_u(i) <= 1 , ...
            alpha_u(i) >= -1, ...
%             alpha_w_2(i) <= 1 , ...
%             alpha_w_2(i) >= -1, ...
%             alpha_w(i) <= 1 , ...
%             alpha_w(i) >= -1,...
            alpha_v(i) <= 1 , ...
            alpha_v(i) >= -1,...
%             alpha_av(i) <= 1 , ...
%             alpha_av(i) >= -1,...            
            ];
    end
    
    %pinv(sys_d.B)*(eye(5)-sys_d.A)*[1;0;0;0;0]

    
    Cost=0;
    for i=1:N
        Cost = Cost + (y{i+1}-ref)'*Qy*(y{i+1}-ref)+ (u{i}-uref)'*Qu*(u{i}-uref);
    end
%    Cost = Cost + (y{i+1}-ref)'*Qy*(y{i+1}-ref);
    options = sdpsettings('verbose',0,'solver','mosek');
    Problem = optimize(Constraints,Cost,options)
    Objective = double(Cost);
    uPred(timesteps) = double(u{1})
    %double(alpha_u(1));
    YPred(:,timesteps) = double(y{2});
    %y11(:,timesteps) = double(y{1});
    norm(YPred(:,timesteps)-ref,2);
    xxx=11;
    
    
 Rplotall{timesteps}= zonotope([ double(R{2}.center), double(R{2}.generators)]);   
    %%  ploting
%     projectedDims = {[1 2],[3 4],[4 5]};
%     for i =1:N+1
%         Rplot{i}= zonotope([ double(R{i}.center), double(R{i}.generators)]) ;
%         Rint{i} = interval(Rplot{i});
%         if i<N+1
%             yplot{i} =double(y{i});
%             uplot{i} =double(u{i});
%         end
%     end
    
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
    %take prev v
    % sep x and y 
    y_t(:,timesteps+1) = sys_d.A * y_t(:,timesteps) + sys_d.B * uPred(timesteps) + w_point +v_point - sys_d.A *prev_v_point;
    yt2ref(timesteps)= norm(y_t(:,timesteps)-ref,2)
    prev_v_point = v_point;
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
            haldy_t= plot(y_t(projectedDims{plotRun}(1),i),y_t(projectedDims{plotRun}(2),i),'*-k');
%         else
%             haldy_t= plot(y0(projectedDims{plotRun}(1)),y0(projectedDims{plotRun}(2)),'*-k');
%         end
    end
    % handC=plot(consSet,projectedDims{plotRun},'k');
    legend([handR,haldy_t],'Reachable set $\mathcal{R}_k$','System trajectory $y(k)$','Interpreter','latex')
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


%sys_d.A*ref+sys_d.B*double(u{1})-ref
%pinv(sys_d.B)*(eye(5)-sys_d.A)*[1;0;0;0;0]

figure
box on;
han_pred=plot(uPred,'b*-');
hold on
han_uref=plot(uref*ones(size(uPred)),'k-');
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
plot(yt2ref,'b*-');
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