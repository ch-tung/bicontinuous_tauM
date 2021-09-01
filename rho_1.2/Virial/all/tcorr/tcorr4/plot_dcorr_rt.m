clear;
close all;

filetype = 'vp';
corrtype = '10';
covtype = 0;

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

%% start loop 
for T_i = ([0.4:0.1:1.0 1.2 1.5 2.0 3.0 5.0]) % loop over temperature
% for T_i = 1.0 % loop over temperature
    tic
    
    iT = iT+1;
    disp(['T = ',num2str(T_i,'%.1f')])
    
    % load dcorr_py.mat file
    filename_dcorr = ['./tcorr4_',filetype,'_py_',num2str(T_i,'%.1f'),'.mat'];
    load(filename_dcorr)
    sim_r_c = eval(['sim_r_',corrtype,'_c']);
    
    % assign color representing the temperature
    cp = parula(20);
    
%     s_fluc = sqrt(sim_s_c(:,12));
    s_fluc = 0*sqrt(1e-2);
    sim_r_normalized = (sim_r_c-sim_a_c)./(sim_s_c + s_fluc^2);
    
    if covtype==1
        sim_r_normalized = (sim_r_c);
    end

    
    % figures
    %% C(r,t)
    figure(iT)
    box on
    log_time = log(ts(2:end));
    [X,Y] = meshgrid(rr(rr>=1),ts(index_t));
    %     sim_r_normalized(sim_r_normalized<1e-8) = 1e-8;
    pcolor(X,Y/1e3,sim_r_normalized(rr>=1,(index_t))')
    %     pcolor(X,Y,(sim_r_c(:,(2:end))'))
    hold on
    %     contour(X,Y,(sim_r_normalized(:,(2:end))'),'w')
%     contour(X,Y,sim_r_normalized(:,(index_t))',exp([0:-1:-6]),'w')
    
%     for n = 0:7
%         v = 10^(-n);
%         plot(0:0.1:10,(0:0.1:10)/v,'--w')
%     end
    
    shading interp
    set(gca, 'YScale', 'log')
    %     set(gca,'ColorScale','log')
    xlim([1 5])
    ylim([1e0 1e7]/1e3)
    xticks([1:10])
    yticks(10.^([-3,0,3]))
    caxis([-0.03 0.06])
    if covtype==1
%         caxis([-0.001 0.001])
    end
    %     colorbar
    
%     colormap(jet)
      
    xlabel('$r$','FontSize',24,'Interpreter','latex')
    ylabel('$t$','FontSize',24,'Interpreter','latex')
    set(gca,'LineWidth',2)
    set(gcf,'Position',[200,100,600,600])
    set(gca,'FontSize',28,'FontName','Times New Roman')
    
    %     saveas(gcf,['C_rt_',num2str(T_i,'%.2f'),'.fig'])
    mkdir(['C_rt_',filetype])
    saveas(gcf,['./C_rt_',filetype,'/C_rt_',corrtype,'_',num2str(T_i,'%.2f'),'.tif'])
       
    toc
end% end of temperature loop
% close all;