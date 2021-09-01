clear
close all

n_particle = 6912;
N = n_particle;       % number of particles
n_conf = 8;
n_grid = 64;

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
% ic = 1;

% index_time = [2 3 6 11 12 15 20 21 24 29 30 33 38 39 42 47 48 51 56 57 60 65];
% index_time = [11];
% index_time = 1:length(t);
% index_time = 21;

%% space
dr = 0.25;
rr = dr:dr:10;

%% temperature
load rho_r_1.2.mat

color = flipud(parula(100));
index_color = round(0.44./Temperature*100);

%% start loop
for T_i = Temperature % loop over temperature
    tic
    
    iT = iT+1;
    disp(['T = ',num2str(T_i,'%.1f')])
    
    %%
    figure(1)
    hold on
    
    Pr_sm = smoothdata(Pr_c(rr>=1,iT),'Gaussian',5);
    plot(rr(rr>=1),Pr_sm/sum(Pr_sm*dr),'Color',color(index_color(iT),:),'LineWidth',2)
%     set(gca, 'XScale', 'log')
%     set(gca, 'YScale', 'log')
    
%     xlabel('\it{r}','FontSize',24)
%     ylabel('\rho({\itr})','FontSize',24)
    
    box on
    
    xlim([0 10])
    xticks([0 5 10])
%     ylim([1e-4 1])
    
    set(gca,'LineWidth',2)
    set(gcf,'Position',[200,100,600,600])
    set(gca,'FontSize',28,'FontName','Arial')
    
    set(gcf,'Position',[200,100,300,300])
    set(gca,'FontSize',24,'FontName','Arial')
    
%     pause(0.1)
    
    Rg(iT) = sqrt(sum(Pr_c(rr<5,iT).*rr(rr<5).^4)/sum(Pr_c(rr<5,iT).*rr(rr<5).^2));
       
    toc
end


