clear
close all

%% load
load('./tau_M.mat')
load('./tcorr_o_6912.mat')
tau_M = tau_M(2:end);
Temperature = T;
Ct = reshape(Ct(:,3,:),[length(t),11]);

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
%     plot(t/tau_M(i),Ct(:,i),'LineWidth',2)
    plot(t/1000,Ct(:,i),'LineWidth',2)
end

% plot([1 1],[0 1],'-k')

set(gca, 'XScale', 'log')

xlim([1e-3 1e4])
xticks(10.^[-3:1:4])
ylim([0 1])

% xlabel('$t/\tau_M$','FontSize',24,'Interpreter','latex')
xlabel('{\itt}','FontSize',24)

ylabel('{\itC}({\itt})','FontSize',24)
set(gca,'LineWidth',2)
set(gcf,'Position',[200,100,600,600])
set(gca,'FontSize',24,'FontName','Arial')
set(gca,'position',[0.1886    0.1874    0.7164    0.7376])

%% susceptibility
load susceptibility_1.2_1.0.mat
for i = 1:length(T)
    [min_dtmax, index_tmax] = min(abs(t-tmax(i)));
    plot(t(index_tmax)/1000,Ct(index_tmax,i),'^','MarkerSize',8,'LineWidth',2,...
        'MarkerFaceColor',color_order(i,:),'MarkerEdgeColor','k')
end

for i = 1:length(T)
    [min_dtau_M, index_tau_M] = min(abs(t-tau_M(i)));
    plot(t(index_tau_M)/1000,Ct(index_tau_M,i),'o','MarkerSize',8,'LineWidth',2,...
        'MarkerFaceColor',color_order(i,:),'MarkerEdgeColor','k')
end

ax = gca;  % grabs current axis
set(ax,'XMinorTick','on')  % sets Minor X Ticks to display