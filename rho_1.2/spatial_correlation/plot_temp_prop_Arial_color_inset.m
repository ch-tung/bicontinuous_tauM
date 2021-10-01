clear
close all

load('curvature_10.mat')
load('tau_M.mat')
tau_M = tau_M(2:end);
T = T(2:end);
v = v(2:end);


%% Figures
% linecolor
color_parula = flipud(parula(100));
index_color = round(0.44./Temperature*100);

% color_order = [0.97 0.58 0.02; 0 0 1];
color_order = [0 0 0; 0 0 0];

figure(1)
hold on
box on
set(gca, 'ColorOrder', color_order)

% yyaxis left

% % plot time chi_max
% load susceptibility_1.2_1.0.mat
% plot(0.44./T,tmax/1000,'^','MarkerSize',12,'LineWidth',2,'Color','#F89406','MarkerFaceColor','#F89406')

%figure settings
box on
set(gca, 'YScale', 'log')
xlim([0 1])
ylim([10 1e5]/1000)

xticks(0:0.2:1)
yticks([1e-2 1e-1 1 1e1 1e2 1e3 1e4])

% xlabel('1/T','FontSize',24)
% ylabel('$\tau_{M}$','FontSize',24,'Interpreter','latex')
set(gca,'LineWidth',2)
set(gcf,'Position',[2000,100,400,400])

%% fit
T_range = T>0.8;

A = 0.44./T(T_range)';
B = log(tau_M(T_range))';

X = [ones(size(A)),A];
Beta = (X'*X)\X'*B;

T_f = 0:0.01:1;
tau_f = exp(Beta(2)*T_f + Beta(1));

% yyaxis left

color_parula = flipud(parula(100));
index_color = round(0.44./T*100);
color_order = color_parula(index_color,:);

plot(T_f,tau_f/1000,'k-')

for i = 1:length(T)
    plot(0.44./T(i),tau_M(i)/1000,'o','MarkerSize',12,'LineWidth',2,...
        'MarkerFaceColor',color_order(i,:),'MarkerEdgeColor','k')
end

%% domain size
% for i = 1:length(Temperature)
%     F_tau = griddedInterpolant(t_MC,mean_GC_T(:,i),'linear');
%     GC_tau_M(i) = F_tau(tau_M(i));
%     k2(i) = -6*GC_tau_M(i);
%     lambda(i) = 2*pi./sqrt(k2(i));
% end
% 
% yyaxis right
% plot(0.44./T,lambda,'sb','MarkerSize',12,'LineWidth',2)
% ylim([0 10])
% ylabel('$\lambda$','FontSize',24,'Interpreter','latex')
% set(gca,'FontSize',28,'FontName','Times New Roman')
% 
%% other densities
% % load lambda_1.3.mat
% % plot(0.44./T,lambda,'sc','MarkerSize',12,'LineWidth',2)
% % ylim([0 10])
% 
% load lambda_1.4.mat
% plot(0.93./Temperature,lambda_GC,'sc','MarkerSize',12,'LineWidth',2)
% ylim([0 10])
% 
% yyaxis left
% plot(0.93./Temperature,tau_M,'o','MarkerSize',12,'LineWidth',2, 'Color', 'r')
% 
% % fit
% T_range = Temperature>2.4;
% 
% A = 0.93./Temperature(T_range)';
% B = log(tau_M(T_range))';
% 
% X = [ones(size(A)),A];
% Beta = (X'*X)\X'*B;
% 
% T_f = 0:0.01:1;
% tau_f = exp(Beta(2)*T_f + Beta(1));
% 
% yyaxis left
% plot(T_f,tau_f,'k-')
%%
% plot(0.44./T,tau_M.*v,'sb-','MarkerSize',12,'LineWidth',2)
% ylim([0 25])

% %% orientation_correlation
% load('./tcorr_o_6912.mat')
% Ct = reshape(Ct(:,3,:),[length(t),11]);
% for i = 1:length(Temperature)
%     F_ct = griddedInterpolant(t,Ct(:,i),'linear');
%     Ct_tau_M(i) = F_ct(tau_M(i));
% end
% 
% yyaxis right
% plot(0.44./T,Ct_tau_M*10,'sg','MarkerSize',12,'LineWidth',2)
% ylim([0 10])
% ylabel('$\lambda$','FontSize',24,'Interpreter','latex')
set(gca,'FontSize',28,'FontName','Times New Roman')

%%
% yyaxis right
% ylim([2 8])
% % ylabel('$\lambda$','FontSize',24,'Interpreter','latex')
% set(gca,'FontSize',28,'FontName','Arial')
% 
% load rho_r_1.2.mat
% plot(0.44./Temperature,Rg,'sb','MarkerSize',12,'LineWidth',2,'MarkerFaceColor','b')

% load rho_r_1.4.mat
% plot(0.93./Temperature,Rg,'xc','MarkerSize',12,'LineWidth',2)
% 
% load rho_r_1.6.mat
% plot(1.76./Temperature,Rg,'^g','MarkerSize',12,'LineWidth',2)
set(gca,'position',[0.1886    0.1874    0.7164    0.7376])

