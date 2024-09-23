clear
close all

load('curvature_10_std.mat')
load('Sigma_interp.mat')

%% Figures_H
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

plot([0 0], [0 1], '--', 'Color', '#666666')
plot([-1 0.5], [0.2 0.2], '--', 'Color', '#666666')

% % fit
% x_fit = mean_MC_T(:,1);
% y_fit = C_c_interp(:,1);
% fo = fitoptions('Method','NonlinearLeastSquares',...
%     'Lower',     [0.1,-0.5],...
%     'Upper',     [10,0.5],...
%     'StartPoint',[1,0]);
% ft = fittype('exp(-exp(a*(x-b)))','options',fo);
% [curve,gof] = fit(x_fit,y_fit,ft);
% Cf = coeffvalues(curve);
% H_list = -1:0.05:1;
% plot(curve)

for i = 1:length(Temperature)
%     patch([mean_MC_T(:,i)-std_MC_T(:,i);flipud(mean_MC_T(:,i)+std_MC_T(:,i))],...
%           [C_c_interp(:,i);flipud(C_c_interp(:,i))],...
%           color_order(i,:),'FaceAlpha',0.1,'LineStyle','none')
    plot(mean_MC_T(:,i),C_c_interp(:,i),'Color',color_order(i,:),'LineWidth',2)
end


ylim([0 1])
xlim([-1 0.5])

xlabel('$\langle H\rangle$','FontSize',28,'Interpreter','latex')
ylabel('$\it\Sigma(t)$','FontSize',24,'Interpreter','latex')

set(gca,'LineWidth',2)
set(gcf,'Position',[200,100,600,600])
set(gca,'FontSize',28,'FontName','Arial')

% ylabel('Averaged mean curvature','FontSize',26)

%% Figures_K
% linecolor
color_parula = flipud(parula(100));
index_color = round(0.44./Temperature*100);

color_order = color_parula(index_color,:);

figure(2)
hold on
box on
set(gca, 'ColorOrder', color_order)

load('tau_M.mat')
tau_M = tau_M(2:end);
tau_T = 1./sqrt(3*Temperature);

for i = 1:length(Temperature)
%     patch([mean_GC_T(:,i)-std_MC_T(:,i);flipud(mean_GC_T(:,i)+std_MC_T(:,i))],...
%           [C_c_interp(:,i);flipud(C_c_interp(:,i))],...
%           color_order(i,:),'FaceAlpha',0.1,'LineStyle','none')
    plot(mean_GC_T(:,i),C_c_interp(:,i),'Color',color_order(i,:),'LineWidth',2)
end


% ylim([0 1])
% xlim([-1 0.5])

xlabel('$\langle K\rangle$','FontSize',28,'Interpreter','latex')
ylabel('$\it\Sigma(t)$','FontSize',24,'Interpreter','latex')

set(gca,'LineWidth',2)
set(gcf,'Position',[200,100,600,600])
set(gca,'FontSize',28,'FontName','Arial')

% ylabel('Averaged mean curvature','FontSize',26)
