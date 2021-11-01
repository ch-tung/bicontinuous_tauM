clear
close all

load('curvature_10.mat')
load('tau_M.mat')
tau_M = tau_M(1:end);
T = T(1:end);
v = v(1:end);


%% Figures
% linecolor
color_parula = flipud(parula(100));
index_color = round(950./T*100);

% color_order = [0.97 0.58 0.02; 0 0 1];
color_order = [0 0 0; 0 0 0];

figure(1)
hold on
box on
set(gca, 'ColorOrder', color_order)

yyaxis left

% % plot time chi_max
% load susceptibility_1.2_1.0.mat
% plot(0.44./T,tmax/1000,'^','MarkerSize',12,'LineWidth',2,'Color','#F89406','MarkerFaceColor','#F89406')

%figure settings
box on
set(gca, 'YScale', 'log')
xlim([0 8e-4])
ylim([1e1 1e3])

xticks(0:2e-4:8e-4)
yticks([1e-2 1e-1 1 1e1 1e2 1e3 1e4])

xlabel('1/T','FontSize',24)
% ylabel('$\tau_{M}$','FontSize',24,'Interpreter','latex')
set(gca,'LineWidth',2)
set(gcf,'Position',[2000,100,600,600])

%% fit
% T_range = T>0.8;
% 
% A = 0.44./T(T_range)';
% B = log(tau_M(T_range))';
% 
% X = [ones(size(A)),A];
% Beta = (X'*X)\X'*B;
% 
% T_f = 0:0.01:1;
% tau_f = exp(Beta(2)*T_f + Beta(1));

yyaxis left
% plot(T_f,tau_f/1000,'k-')

plot(1./T,tau_M,'o','MarkerSize',12,'LineWidth',2,'Color','r','MarkerFaceColor','r')

%% domain size
% for i = 1:length(T)
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
% plot(0.93./T,lambda_GC,'sc','MarkerSize',12,'LineWidth',2)
% ylim([0 10])
% 
% yyaxis left
% plot(0.93./T,tau_M,'o','MarkerSize',12,'LineWidth',2, 'Color', 'r')
% 
% % fit
% T_range = T>2.4;
% 
% A = 0.93./T(T_range)';
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
% for i = 1:length(T)
%     F_ct = griddedInterpolant(t,Ct(:,i),'linear');
%     Ct_tau_M(i) = F_ct(tau_M(i));
% end
% 
% yyaxis right
% plot(0.44./T,Ct_tau_M*10,'sg','MarkerSize',12,'LineWidth',2)
% ylim([0 10])
% ylabel('$\lambda$','FontSize',24,'Interpreter','latex')
% set(gca,'FontSize',28,'FontName','Times New Roman')

%%
yyaxis right
ylim([0 10])
% ylabel('$\lambda$','FontSize',24,'Interpreter','latex')
set(gca,'FontSize',28,'FontName','Arial')

load rho_r_Fe.mat
plot(1./T,Rg,'sb','MarkerSize',12,'LineWidth',2,'MarkerFaceColor','b')

% load rho_r_1.4.mat
% plot(0.93./T,Rg,'xc','MarkerSize',12,'LineWidth',2)
% 
% load rho_r_1.6.mat
% plot(1.76./T,Rg,'^g','MarkerSize',12,'LineWidth',2)
set(gca,'position',[0.1886    0.1874    0.7164    0.7376])

