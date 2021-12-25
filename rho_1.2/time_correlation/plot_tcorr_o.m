clear
close all

%% load
load('./tau_M.mat')
load('./tcorr_o_6912_extreme.mat')
tau_M = tau_M(2:end);
Temperature = T;
Ct = reshape(Ct(:,3,:),[length(t),11]);

%% Figures
% linecolor
color_parula = flipud(parula(100));
index_color = round(0.44./Temperature*100);

color_order = color_parula(index_color,:);

figure(1)
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
set(gca,'FontSize',28,'FontName','Arial')
set(gca,'position',[0.1886    0.1874    0.7376    0.7376])

% susceptibility
load susceptibility_1.2_1.0.mat
for i = 1:length(T)
    F_tau = griddedInterpolant(t,Ct(:,i),'linear');
    Ct_tau_M(i) = F_tau(tmax(i));
%     [min_dtmax, index_tmax] = min(abs(t-tmax(i)));
    plot(tmax(i)/1000,Ct_tau_M(i),'^','MarkerSize',8,'LineWidth',2,...
        'MarkerFaceColor',color_order(i,:),'MarkerEdgeColor','k')
end

for i = 1:length(T)
    F_tau = griddedInterpolant(t,Ct(:,i),'linear');
    Ct_tau_M(i) = F_tau(tau_M(i));
    [min_dtau_M, index_tau_M] = min(abs(t-tau_M(i)));
    plot(tau_M(i)/1000,Ct_tau_M(i),'o','MarkerSize',8,'LineWidth',2,...
        'MarkerFaceColor',color_order(i,:),'MarkerEdgeColor','k')
end

ax = gca;  % grabs current axis
set(ax,'XMinorTick','on')  % sets Minor X Ticks to display

%%
figure(2)
hold on
box on
set(gca, 'ColorOrder', color_order)

for i = 1:length(T)
%     plot(t/tau_M(i),Ct(:,i),'LineWidth',2)
    plot(t/1000,Ct_1(:,3,i),'LineWidth',2)
end
for i = 1:length(T)
%     plot(t/tau_M(i),Ct(:,i),'LineWidth',2)
    plot(t/1000,Ct_0(:,3,i),'LineWidth',2)
end


% plot([1 1],[0 1],'-k')

set(gca, 'XScale', 'log')

xlim([1e-3 1e4])
xticks(10.^[-3:1:4])
% ylim([0 1])

% xlabel('$t/\tau_M$','FontSize',24,'Interpreter','latex')
xlabel('{\itt}','FontSize',24)

ylabel('{\itC}({\itt})','FontSize',24)
set(gca,'LineWidth',2)
set(gcf,'Position',[200,100,600,600])
set(gca,'FontSize',28,'FontName','Arial')
set(gca,'position',[0.1886    0.1874    0.7376    0.7376])