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

color_order = [1 0 0; 0 0 1];

figure(1)
hold on
box on
set(gca, 'ColorOrder', color_order)

yyaxis left
plot(0.44./T,tau_M/1000,'or','MarkerSize',12,'LineWidth',2)

%figure settings
box on
set(gca, 'YScale', 'log')
xlim([0 1])
ylim([10 1e5]/1000)

xticks(0:0.2:1)

xlabel('$T_g/T$','FontSize',24,'Interpreter','latex')
ylabel('$\tau_{M}$','FontSize',24,'Interpreter','latex')
set(gca,'LineWidth',2)
set(gcf,'Position',[2000,100,800,600])

%% fit
T_range = T>0.8;

A = 0.44./T(T_range)';
B = log(tau_M(T_range))';

X = [ones(size(A)),A];
Beta = (X'*X)\X'*B;

T_f = 0:0.01:1;
tau_f = exp(Beta(2)*T_f + Beta(1));

yyaxis left
plot(T_f,tau_f/1000,'k-')

%% domain size
for i = 1:length(Temperature)
    F_tau = griddedInterpolant(t_MC,mean_GC_T(:,i),'linear');
    GC_tau_M(i) = F_tau(tau_M(i));
    k2(i) = -6*GC_tau_M(i);
    lambda(i) = 2*pi./sqrt(k2(i));
end

yyaxis right
plot(0.44./T,lambda,'xb','MarkerSize',12,'LineWidth',2)
ylim([0 10])
ylabel('$\lambda$','FontSize',24,'Interpreter','latex')
set(gca,'FontSize',28,'FontName','Times New Roman')

% plot(0.44./T,tau_M.*v,'sb-','MarkerSize',12,'LineWidth',2)
% ylim([0 25])

%% orientation_correlation
load('./tcorr_o_6912.mat')
Ct = reshape(Ct(:,3,:),[length(t),11]);
for i = 1:length(Temperature)
    F_ct = griddedInterpolant(t,Ct(:,i),'linear');
    Ct_tau_M(i) = F_ct(tau_M(i));
end

yyaxis right
plot(0.44./T,Ct_tau_M*10,'sg','MarkerSize',12,'LineWidth',2)
ylim([0 10])
ylabel('$\lambda$','FontSize',24,'Interpreter','latex')
set(gca,'FontSize',28,'FontName','Times New Roman')