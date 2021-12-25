clear;
close all;

n_conf = 8;
N = 6912;

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
% index_time = [2 3 6 11 12 15 20 21 24 29 30 33 38 39 42 47 48 51 56 57 60 65];
index_time = 1:length(t);
%% start loop
T = [0.5:0.1:1.0 1.2 1.5 2.0 3.0 5.0];

color_parula = flipud(parula(100));
index_color = round(0.44./T*100);
color_order = color_parula(index_color,:);

for T_i = T % loop temperature
    
    iT = iT+1;
    disp(['T = ',num2str(T_i,'%.1f')])
    filename_rg = ['../rg_mat/rg_py_',num2str(T_i,'%.1f'),'.mat'];
    load(filename_rg)
    
    %% t = 0
    for it = [1] % loop time
        coord_0 = [];
        p_0 = [];
        Vs_0_C = [];
        for ic = 1:n_conf % loop conf
            %                 disp(['C = ',num2str(ic,'%.0f')])
            rg_all = rg_all_t_c(:,:,:,it,ic);
            [p,coord] = calculate_shape(rg_all);
            p_0 = [p_0;p];
            coord_0 = [coord_0;coord];
            for i1 = 1:size(rg_all,3)
                [V,D] = eig(rg_all(:,:,i1));
                [lambda,ind] = sort(diag(D));
                Vs(:,:,i1) = V(:,ind);
            end
            Vs_0_C = cat(3,Vs_0_C,Vs);
        end
    end
    
    %% loop over t
    for it = index_time % loop
        Vs_C = [];
        disp(['t = ',num2str(t(it))])
        for ic = 1:n_conf % loop
            %                 disp(['C = ',num2str(ic,'%.0f')])
            rg_all = rg_all_t_c(:,:,:,it,ic);
            for i1 = 1:size(rg_all,3)
                [V,D] = eig(rg_all(:,:,i1));
                [lambda,ind] = sort(diag(D));
                Vs(:,:,i1) = V(:,ind);
            end
            Vs_C = cat(3,Vs_C,Vs);
        end
               
        sl0 = 1:N*n_conf;
        
        n_extreme = round(N/10);
        for j = 1:3
            ct = (3*(dot(Vs_0_C(:,j,sl0),Vs_C(:,j,sl0))).^2-1)/2;
            Ct(it,j,iT) = mean(ct,3);
            Ct_var(it,j,iT) = var(ct);
            
            [ct_sort,i_sort] = sort(ct(:));
            Ct_1(it,j,iT) = mean(ct_sort(end-n_extreme+1:end));
            Ct_0(it,j,iT) = mean(ct_sort(1:n_extreme));
        end
        
    end
    %         Ct = Ct./Ct(1,:);
    figure(1)
    hold on
    plot(t,Ct(:,3,iT),'-','LineWidth',2)
    plot(t,Ct(:,3,iT)-sqrt(Ct_var(:,3,iT)),':','LineWidth',2)
    plot(t,Ct(:,3,iT)+sqrt(Ct_var(:,3,iT)),':','LineWidth',2)
    plot(t,Ct_0(:,3,iT),'--','LineWidth',1)
    plot(t,Ct_1(:,3,iT),'--','LineWidth',1)
    set(gca, 'XScale', 'log')
    
end

%% figure settings
box on
set(gca,'LineWidth',2)
xlabel('$t$','FontSize',24,'Interpreter','latex')
ylabel('$C(t)$','FontSize',24,'Interpreter','latex')

% set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
set(gca, 'LineWidth', 2)
xlim([1e2 1e7])
ylim([0 1])

set(gca,'FontSize',28,'FontName','Times New Roman')
set(gcf,'Position',[2000,100,800,600])