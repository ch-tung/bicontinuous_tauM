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
tau_M = tau_M(2:end)/1000;

t_MC = t_MC/1000;

plot(t_MC,mean_MC_T,'LineWidth',2)

for i = 1:length(Temperature)
    F_tau = griddedInterpolant(t_MC,mean_MC_T(:,i),'linear');
    MC_tau_M(i) = F_tau(tau_M(i));
    plot(tau_M(i),MC_tau_M(i),'o','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor',color_order(i,:),...
        'MarkerSize',8)
end

load susceptibility_1.2_1.0.mat
for i = 1:length(Temperature)
    F_max = griddedInterpolant(t_MC,mean_MC_T(:,i),'linear');
    MC_tmax(i) = F_max(tmax(i)/1000);
    plot(tmax(i)/1000,MC_tmax(i),'^','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor',color_order(i,:),...
        'MarkerSize',8)
end

plot([1e-5 1e5],[0 0],'-','Color','#666666')

set(gca, 'XScale', 'log')
xlim([1e-3 1e4])
ylim([-0.5 0.5])
xticks(10.^[-3:1:4])

xlabel('\it{t}','FontSize',28)

set(gca,'LineWidth',2)
set(gcf,'Position',[200,100,600,600])
set(gca,'FontSize',24,'FontName','Arial')

ylabel('<H>','FontSize',24)
set(gca,'position',[0.1886    0.1874    0.7164    0.7376])

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
% xlim([1e-3 1e4])
% ylim([0 20])
% xticks([1e-3 1 1e3])
% 
% xlabel('$t$','FontSize',24,'Interpreter','latex')
% ylabel('$2\pi\left\langle k^2 \right\rangle^{-\frac{1}{2}}$','FontSize',24,'Interpreter','latex')
% set(gca,'LineWidth',2)
% set(gcf,'Position',[200,100,800,600])
% set(gca,'FontSize',28,'FontName','Times New Roman')
