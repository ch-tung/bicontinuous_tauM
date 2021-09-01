clear
close all

clear;
close all;

filetype = 'vs';

%% time
dt = 1e-2;
q = 7.25;
n_ave = 1;
n_inteval = 1;
n_block = n_ave*n_inteval;

ts = [0 1:9];
for n_2 = 1:7
    ts = [ts n_ave*(10^n_2)*[1:9]];
end
ts = ts*n_inteval;
ts = ts(ts<=1e7);

n_frame = length(ts);

t = ts;

iT = 0;
index_t = [2 3 6 11 12 15 20 21 24 29 30 33 38 39 42 47 48 51 56 57 60 65];
% index_t = 2:2:length(t);
% index_t = 1:length(t);

%% Temperature
Temperature = fliplr([0.5:0.1:1.0 1.2 1.5 2.0 3.0 5.0]);

% linecolor
color = flipud(parula(100));
index_color = round(0.44./Temperature*100);

% load('YlOrRd.mat')
% color = rgba(:,1:3);
% index_color = round((0.44./Temperature*0.75+0.25)*length(color));

    
load tau_M.mat
tau_M = tau_M(T>=0.5);
tau_M = tau_M/1000;
tau_M = fliplr(tau_M);

%% start loop
for T_i = Temperature % loop over temperature
    % for T_i = 1.0 % loop over temperature
    tic
    
    iT = iT+1;
    disp(['T = ',num2str(T_i,'%.1f')])
    
    % load dcorr_py.mat file
    filename_dcorr = ['./tSD_CG_ex_',filetype,'_py_',num2str(T_i,'%.1f'),'.mat'];
    load(filename_dcorr)
    
    %%
    figure(1)
    hold on
    box on
    
    t = ts/1000;
    plot([1e-6 1e6],[0 0], '-','Color','#666666')
    
    C_0 = smoothdata(C_0,'gaussian',5);
    C_1 = smoothdata(C_1,'gaussian',5);
    C = smoothdata(C,'gaussian',5);
    
%     plot(t(index_t),(C_0(index_t)-C(index_t))/mean(C(60:end)),'-','Color',color(index_color(iT),:),'LineWidth',2)
%     plot(t(index_t),(C_1(index_t)-C(index_t))/mean(C(60:end)),'--','Color',color(index_color(iT),:),'LineWidth',2)
    
    plot(t(index_t),C_0(index_t)/mean(C(60:end)),'-','Color',color(index_color(iT),:),'LineWidth',2)
    plot(t(index_t),C_1(index_t)/mean(C(60:end)),'--','Color',color(index_color(iT),:),'LineWidth',2)
    
    set(gca, 'XScale', 'log')
    
    xlim([1 1e7]/1000)
    xticks(10.^[-3:1:4])
    ylim([0 1.2])
    yticks(0:0.2:1.2)
    
%     yrange = 0.04;
%     ylim([-yrange yrange])
%     yticks([-yrange 0 yrange])
    
    xlabel('\it{t}','FontSize',24)
%     ylabel('$SD-\left\langle SD \right\rangle$','FontSize',24,'Interpreter','latex')
    ylabel('\it{SD}','FontSize',24)
    set(gca,'LineWidth',2)
    set(gcf,'Position',[200,100,600,600])
    set(gca,'FontSize',28,'FontName','Arial')
    set(gca,'position',[0.1886    0.1874    0.7164    0.7376])
    
    ax=gca;
    ax.XAxis.FontSize = 24;
    
    %     saveas(gcf,['C_rt_',num2str(T_i,'%.2f'),'.fig'])
    %     mkdir(['C_rt_',filetype])
    %     saveas(gcf,['./C_rt_',filetype,'/C_rt_',num2str(T_i,'%.2f'),'.tif'])
    
    %%
    figure(2)
    hold on
    box on
    
    t = ts/1000;
    plot([1e-6 1e6],[0 0], '-','Color','#666666')
    plot([1 1],[-1 1], '--','Color','#666666')
    
    C_0 = smoothdata(C_0,'gaussian',5);
    C_1 = smoothdata(C_1,'gaussian',5);
    C = smoothdata(C,'gaussian',5);
    
    plot(t(index_t)/tau_M(iT),(C_0(index_t)-C(index_t))/mean(C(60:end)),'-','Color',color(index_color(iT),:),'LineWidth',2)
    plot(t(index_t)/tau_M(iT),(C_1(index_t)-C(index_t))/mean(C(60:end)),'--','Color',color(index_color(iT),:),'LineWidth',2)
    
%     plot(t(index_t),C_0(index_t)/mean(C(60:end)),'-','Color',color(index_color(iT),:),'LineWidth',2)
%     plot(t(index_t),C_1(index_t)/mean(C(60:end)),'--','Color',color(index_color(iT),:),'LineWidth',2)
    
    set(gca, 'XScale', 'log')
    
    xlim(10.^[-5 5])
    xticks(10.^[-10:2:10])
    
    yrange = 0.02;
    ylim([-yrange yrange])
    yticks([-yrange 0 yrange])
    
    xlabel('\it{t}/\it{\tau_M}','FontSize',24)
    ylabel('$SD-\left\langle SD \right\rangle$','FontSize',24,'Interpreter','latex')
%     ylabel('\it{SD}','FontSize',24)
    set(gca,'LineWidth',2)
    set(gcf,'Position',[200,100,600,600])
    set(gca,'FontSize',28,'FontName','Arial')
    set(gca,'position',[0.1886    0.1874    0.7164    0.7376])
    
    %     saveas(gcf,['C_rt_',num2str(T_i,'%.2f'),'.fig'])
    %     mkdir(['C_rt_',filetype])
    %     saveas(gcf,['./C_rt_',filetype,'/C_rt_',num2str(T_i,'%.2f'),'.tif'])
    
    toc
end% end of temperature loop
% close all;
