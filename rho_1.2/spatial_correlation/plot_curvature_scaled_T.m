clear
close all

load('curvature_10.mat')

%% Figures
% linecolor
color_parula = flipud(parula(100));
index_color = round(0.44./Temperature*100);

color_order = color_parula(index_color,:);

figure(1)
hold on
box on
set(gca, 'ColorOrder', color_order)

load('tau_M.mat')
tau_M = tau_M(2:end);
tau_T = 1./sqrt(3*Temperature);

plot(t_MC./tau_T/1e3,mean_MC_T,'LineWidth',2)

% for i = 1:length(Temperature)
%     F_tau = griddedInterpolant(t_MC,mean_MC_T(:,i),'linear');
%     MC_tau_M(i) = F_tau(tau_M(i));
%     plot(1,MC_tau_M(i),'o','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor',color_order(i,:),...
%         'MarkerSize',12)
% end

plot([1e-5 1e5],[0 0],'-','Color','#666666')
plot([1 1],[-1 1],'--','Color','#666666')

set(gca, 'XScale', 'log')
xlim(10.^[-4 4])
xticks(10.^[-10:2:10])
ylim([-0.5 0.5])

xlabel('\it{t}/\tau_{M}','FontSize',28)

set(gca,'LineWidth',2)
set(gcf,'Position',[200,100,600,600])
set(gca,'FontSize',28,'FontName','Arial')

% ylabel('Averaged mean curvature','FontSize',26)

%%
figure(1)
hold on
box on

yyaxis right
set(gca, 'ColorOrder', color_order)

load('tau_M.mat')
tau_M = tau_M(2:end);
tau_T = 1./sqrt(Temperature);

plot(t_MC./tau_T/1e3,mean_GC_T,'--','LineWidth',2)

% for i = 1:length(Temperature)
%     F_tau = griddedInterpolant(t_MC,mean_MC_T(:,i),'linear');
%     MC_tau_M(i) = F_tau(tau_M(i));
%     plot(1,MC_tau_M(i),'o','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor',color_order(i,:),...
%         'MarkerSize',12)
% end

plot([1e-8 1e8],[0 0],'-','Color','#666666')
plot([1 1],[-1 1],'--','Color','#666666')

set(gca, 'XScale', 'log')
xlim(10.^[-4 4])
xticks(10.^[-10:2:10])
ylim([-0.5 0.5])

xlabel('\it{t}v_{RMS}','FontSize',28)

set(gca,'LineWidth',2)
set(gcf,'Position',[200,100,600,600])
set(gca,'FontSize',28,'FontName','Arial')

% ylabel('Averaged mean curvature','FontSize',26)
set(gca,'position',[0.13    0.22   0.7376    0.7376])
%%
% figure(2)
% hold on
% box on
% set(gca, 'ColorOrder', color_order)
% 
% mean_GC_T(mean_GC_T>0) = 0;
% 
% plot(t_MC,2*pi./sqrt(-6*mean_GC_T),'LineWidth',2)
% set(gca, 'XScale', 'log')
% 
% for i = 1:length(Temperature)
%     F_tau = griddedInterpolant(t_MC,mean_GC_T(:,i),'linear');
%     GC_tau_M(i) = F_tau(tau_M(i));
%     plot(tau_M(i),2*pi./sqrt(-6*GC_tau_M(i)),'o','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor',color_order(i,:),...
%         'MarkerSize',12)
% end
% 
% xlim([10 1e6])
% ylim([0 20])
% 
% xlabel('$t$','FontSize',24,'Interpreter','latex')
% ylabel('$2\pi\left\langle k^2 \right\rangle^{-\frac{1}{2}}$','FontSize',24,'Interpreter','latex')
% set(gca,'LineWidth',2)
% set(gcf,'Position',[200,100,800,600])
% set(gca,'FontSize',28,'FontName','Times New Roman')
