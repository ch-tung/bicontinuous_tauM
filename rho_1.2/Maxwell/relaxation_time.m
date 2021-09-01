% Calculate Maxwell relaxation time
clear
% close all

n_conf = 32;
% T = [0.1 0.2 0.3];
% T = [0.4 0.5 0.6];
T = [0.4:0.1:1.0 1.2 1.5 2.0 3.0 5.0];
Tc = 0.44;

% linecolor
color_parula = flipud(parula(100));
index_color = round(Tc./T(2:end)*100);
color_order = color_parula(index_color,:);

rho = 1.2;
N = 10976;

filename = 'maxwell_long_py.mat';
load(filename)

t = double(t'-t(1));

%% deal with noise
Cmin = 0.00001;
for i = 1:size(C_c,2)
   index = find(C_c(:,i)<Cmin);
   if isempty(index)
       continue
   end
   i_min(i) = min(index);
   C_c(i_min(i):end,i) = Cmin;
end

figure(1)
hold on
box on
set(gca, 'ColorOrder', color_order)
plot(t,C_c(:,2:end),'LineWidth',2)
set(gca, 'XScale', 'log')

xlim([10 1e6])
ylim([0 1])

xlabel('$t$','FontSize',24,'Interpreter','latex')
set(gca,'LineWidth',2)
set(gcf,'Position',[2000,100,600,600])
set(gca,'FontSize',28,'FontName','Times New Roman')
ylabel('$\left\langle \sigma^{\textrm{xy}}(0)\sigma^{\textrm{xy}}(t) \right\rangle / \left\langle (\sigma^{\textrm{xy}}(0))^{2} \right\rangle$','FontSize',24,'Interpreter','latex')


%% itegral
tau_M = trapz(t,C_c);

figure(2)
hold on

set(gca, 'ColorOrder', color_order)

for i = 2:length(T)
plot(Tc./T(i),tau_M(i),'or','LineWidth',2,'MarkerSize',10)
end

%figure settings
box on
set(gca, 'YScale', 'log')
xlim([0 1])
ylim([40 40000])

xticks(0:0.2:1)

xlabel('$T_g/T$','FontSize',24,'Interpreter','latex')
ylabel('$\tau_{M}$','FontSize',24,'Interpreter','latex')
set(gca,'LineWidth',2)
set(gcf,'Position',[2000,100,600,600])
set(gca,'FontSize',28,'FontName','Times New Roman')

%% fit
T_range = T>0.8;

A = Tc./T(T_range)';
B = log(tau_M(T_range))';

X = [ones(size(A)),A];
Beta = (X'*X)\X'*B;

T_f = 0:0.01:1; % T_f = T_c/T
tau_f = exp(Beta(2)*T_f + Beta(1));

plot(T_f,tau_f,'k-')

slope = -Beta(2)*Tc;
Ea = -slope;
tau_inf = exp(Beta(2)*0 + Beta(1));

%% velocity
c_T = sqrt(G_V_c/(1.2*N));
