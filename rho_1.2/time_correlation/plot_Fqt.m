clear
close all

%% load
load('./tau_M.mat')
load('./Fqt_864.mat')
tau_M = tau_M(2:end);
Temperature = T;
Fqt = Fqtm_s./Fq0m_s;

%% Figures
% linecolor
color_parula = flipud(parula(100));
index_color = round(0.44./Temperature*100);

color_order = color_parula(index_color,:);

figure(2)
hold on
box on
set(gca, 'ColorOrder', color_order)

for i = 1:length(T)
    plot(t/tau_M(i),Fqt(:,i),'LineWidth',2)
%     plot(t/1000,Fqt(:,i),'LineWidth',2)
end

plot([1 1],[0 1],'-k')

set(gca, 'XScale', 'log')

% xlim([1e-3 1e4])
ylim([0 1])

% xlabel('$t/\tau_M$','FontSize',24,'Interpreter','latex')
xlabel('$t$','FontSize',24,'Interpreter','latex')

ylabel('$F_{s}(q,t)$','FontSize',24,'Interpreter','latex')
set(gca,'LineWidth',2)
set(gcf,'Position',[200,100,800,600])
set(gca,'FontSize',28,'FontName','Times New Roman')
