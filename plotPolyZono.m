clear all
close all

%load work spaces
load('workspaces\poly')
load('workspaces\ZPC')

%less noise
%initpoints 4 step 100
%load('workspaces\ZPC_in4st100W0.01V0.002.mat');
%load('workspaces\poly_W0.01V0.002.mat');


%More noise
%load('workspaces\ZPC_in100st5W0.1V0.02N2.mat')
%load('workspaces\poly_W0.1V0.02N2.mat')


%% plot 3d reachable sets over steps
projectedDims3d = {[1 2 3],[3 4 5]};
for plotRun=1:length(projectedDims3d)
    figure('Renderer', 'painters', 'Position', [10 10 900 900]);
    hold on;
    box on;
    for i=1:maxsteps
        handR=  plot(Rplotall{i},projectedDims3d{plotRun},'r','Filled',false);
        % plot(YPred(projectedDims{plotRun}(1),i),YPred(projectedDims{plotRun}(2),i),'*-k');
    end
    for i=2:maxsteps+1
        haldy_t= plot3(y_t(projectedDims3d{plotRun}(1),i),y_t(projectedDims3d{plotRun}(2),i),y_t(projectedDims3d{plotRun}(3),i),'+b');
        %         else
        haldy_Pred= plot3(YPred(projectedDims3d{plotRun}(1),i),YPred(projectedDims3d{plotRun}(2),i),YPred(projectedDims3d{plotRun}(3),i),'*k');
        %         end
    end
    % handC=plot(consSet,projectedDims{plotRun},'k');
    legend([handR,haldy_t,haldy_Pred],'Reachable set $\hat{\mathcal{R}}_t$','System trajectory $y(t)$','$y$-pred$(t)$','Interpreter','latex')
    % label plot
    %xlabel(['x_{',num2str(projectedDims{plotRun}(1)),'}']);
    %ylabel(['x_{',num2str(projectedDims{plotRun}(2)),'}']);
    xlabel(['$y_{',num2str(projectedDims3d{plotRun}(1)),'}$'],'Interpreter','latex');
    ylabel(['$y_{',num2str(projectedDims3d{plotRun}(2)),'}$'],'Interpreter','latex');
    zlabel(['$y_{',num2str(projectedDims3d{plotRun}(3)),'}$'],'Interpreter','latex');
    
    ax = gca;
    ax.FontSize = 19;
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

%% plot 2d reachable sets over steps
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
    legend([handR,haldy_t,haldy_Pred],'Reachable set $\hat{\mathcal{R}}_t$','System trajectory $y(t)$','$y$-pred$(t)$','Interpreter','latex')
    % label plot
    %xlabel(['x_{',num2str(projectedDims{plotRun}(1)),'}']);
    %ylabel(['x_{',num2str(projectedDims{plotRun}(2)),'}']);
    xlabel(['$y_{',num2str(projectedDims{plotRun}(1)),'}$'],'Interpreter','latex');
    ylabel(['$y_{',num2str(projectedDims{plotRun}(2)),'}$'],'Interpreter','latex');
    
    
    ax = gca;
    ax.FontSize = 19;
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

figure('Renderer', 'painters', 'Position', [10 10 700 800]);
box on;
hold on
han_pred=plot(uPred,'b*-');
%using polytopes
han_pred_poly=plot(uk,'r-');
%using zonotopes
han_pred_model=plot(uPred_model,'r+-');

han_uref=plot(uref*ones(size(uPred)),'k-');
ax = gca;
xlabel('Time step $t$','Interpreter','latex')
legend([han_pred,han_pred_model,han_pred_poly,han_uref],'ZPC $u^*(t)$','RMPC-zono $u^*(t)$','RMPC-poly $u^*(t)$','Reference input $r_u(t)$','Interpreter','latex')
ax.FontSize = 19;
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
    legend([handR,haldy_t],'Reachable set $\hat{\mathcal{R}}_t$','$y$-pred$(t)$','Interpreter','latex')
    % label plot
    %xlabel(['x_{',num2str(projectedDims{plotRun}(1)),'}']);
    %ylabel(['x_{',num2str(projectedDims{plotRun}(2)),'}']);
    xlabel(['$y_{',num2str(projectedDims{plotRun}(1)),'}$'],'Interpreter','latex');
    ylabel(['$y_{',num2str(projectedDims{plotRun}(2)),'}$'],'Interpreter','latex');
    
    
    ax = gca;
    ax.FontSize = 19;
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
%% ---------------------- plot y's ZPC----------------------------------------%%
figure('Renderer', 'painters', 'Position', [10 10 700 800]);
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
%ylim([-2.7,5.3])
xlim([0,80])
xlabel('Time step $t$','Interpreter','latex')
legend([handy{1},handy_pred{1},handy{2},handy_pred{2},handy{3},handy_pred{3},handy{4},handy_pred{4},handy{5},handy_pred{5}],...
    '$y_1(t)$','$y_1$-pred$(t)$','$y_2(t)$','$y_2$-pred$(t)$','$y_3(t)$','$y_3$-pred$(t)$','$y_4(t)$','$y_4$-pred$(t)$','$y_5(t)$','$y_5$-pred$(t)$','Interpreter','latex');
ax = gca;
ax.FontSize = 19;
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

%% ---------------------- plot y's MPC----------------------------------------%%
% figure
% hold on
% box on;
% %for i =1:dim_x
% handy{1}= plot(y_t_model(1,:)','r');
% handy_pred{1}= plot(YPred_model(1,:)','-*r');
% handy{2}= plot(y_t_model(2,:)','k');
% handy_pred{2}= plot(YPred_model(2,:)','-*k');
% handy{3}= plot(y_t_model(3,:)','b');
% handy_pred{3}= plot(YPred_model(3,:)','-*b');
% handy{4}= plot(y_t_model(4,:)','g');
% handy_pred{4}= plot(YPred_model(4,:)','-*g');
% handy{5}= plot(y_t_model(5,:)','m');
% handy_pred{5}= plot(YPred_model(5,:)','m-*');
% %end
% 
% % for i =1:dim_x
% %  handy_pred{i}= plot(YPred(i,:)','-*');
% % end
% xlabel('Time step $t$','Interpreter','latex')
% legend([handy{1},handy_pred{1},handy{2},handy_pred{2},handy{3},handy_pred{3},handy{4},handy_pred{4},handy{5},handy_pred{5}],...
%     '$y_1(t)$','$y_1$-pred$(t)$','$y_2(t)$','$y_2$-pred$(t)$','$y_3(t)$','$y_3$-pred$(t)$','$y_4(t)$','$y_4$-pred$(t)$','$y_5(t)$','$y_5$-pred$(t)$','Interpreter','latex');
% ax = gca;
% ax.FontSize = 12;
% %set(gcf, 'Position',  [50, 50, 800, 400])
% ax = gca;
% outerpos = ax.OuterPosition;
% ti = ax.TightInset;
% left = outerpos(1) + ti(1);
% bottom = outerpos(2) + ti(2);
% ax_width = outerpos(3) - ti(1) - ti(3)-0.01;
% ax_height = outerpos(4) - ti(2) - ti(4);
% ax.Position = [left bottom ax_width ax_height];

%legend([handy{1},handy{2},handy{3},handy{4},handy{5}],...
%        'y1','y2','y3','y4','y5','Location','northwest');

%% -------------------------ZPC yt2ref-------------------------------------
figure('Renderer', 'painters', 'Position', [10 10 700 800]);
hold on
box on;
han_yt2ref_poly=plot(yt2ref_poly,'r-');
han_yt2ref_model=plot(yt2ref_model,'r+-');
han_yt2ref=plot(yt2ref,'b*-');
%ylim([0,2.5])
xlabel('Time step $t$','Interpreter','latex')
legend([han_yt2ref,han_yt2ref_model,han_yt2ref_poly],'ZPC $|| y(t) - r_y(t) ||$','RMPC-zono $|| y(t) - r_y(t) ||$','RMPC-poly $|| y(t) - r_y(t) ||$','Interpreter','latex')
ax = gca;
ax.FontSize = 19;
%set(gcf, 'Position',  [50, 50, 800, 400])
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3)-0.01;
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];


%% ------------------------- cost -------------------------------------
figure('Renderer', 'painters', 'Position', [10 10 700 800]);
hold on
box on;
han_yt2ref_poly=plot(Cost_rob_ol,'r-');
han_yt2ref_model=plot(Cost_model_vec,'r+-');
han_yt2ref=plot(Cost_vec,'b*-');
%ylim([0,3000])
xlabel('Time step $t$','Interpreter','latex')
legend([han_yt2ref,han_yt2ref_model,han_yt2ref_poly],'ZPC cost','RMPC-zono cost','RMPC-poly cost','Interpreter','latex')
ax = gca;
ax.FontSize = 19;
%set(gcf, 'Position',  [50, 50, 800, 400])
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3)-0.01;
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

%% %%%%%%%%%%%%%%%%%% constrain on y2%%%%%%%%%%%%%%%%%%%%%%%%



figure('Renderer', 'painters', 'Position', [10 10 700 800]);
hold on
box on;
hand_y_t= plot(y_t(2,:),'b*-');
hand_y_poly= plot(Y(2,:),'r-');
hand_y_t_model= plot(y_t_model(2,:),'r+-');
hand_YPred_t=plot(YPred(2,:),'b+-');

handcon = plot(intc.inf(2)*ones(size(y_t(2,:))),'k--');
handcon = plot(intc.sup(2)*ones(size(y_t(2,:))),'k--');
axis([0,40 ,min(min(y_t(2,:)),intc.sup(2))-1, intc.sup(2)+1]);
xlabel('Time step $t$','Interpreter','latex')
legend([hand_y_t,hand_YPred_t,hand_y_t_model,hand_y_poly,handcon],'ZPC $y_2(t)$','ZPC $y_2$-pred$(t)$','RMPC-zono $y_2(t)$','RMPC-poly $y_2(t)$','constraint','Interpreter','latex');
ax = gca;
ax.FontSize = 19;
%set(gcf, 'Position',  [50, 50, 800, 400])
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3)-0.01;
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
