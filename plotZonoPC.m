clear all
close all
load('workspaces/run100iterworkspace.mat')



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
    xlabel(['$y_{',num2str(projectedDims{plotRun}(1)),'}$'],'Interpreter','latex');
    ylabel(['$y_{',num2str(projectedDims{plotRun}(2)),'}$'],'Interpreter','latex');
    
    
    ax = gca;
    ax.FontSize = 18;
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
    xlabel(['$y_{',num2str(projectedDims{plotRun}(1)),'}$'],'Interpreter','latex');
    ylabel(['$y_{',num2str(projectedDims{plotRun}(2)),'}$'],'Interpreter','latex');
    
    
    ax = gca;
    ax.FontSize = 18;
    %set(gcf, 'Position',  [50, 50, 800, 400])
    ax = gca;
    outerpos = ax.OuterPosition;
    ti = ax.TightInset;
    left = outerpos(1) + ti(1);
    bottom = outerpos(2) + ti(2);
    ax_width = outerpos(3) - ti(1) - ti(3)-0.02;
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