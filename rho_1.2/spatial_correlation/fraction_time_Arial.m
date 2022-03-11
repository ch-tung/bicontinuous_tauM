clear
close all

n_particle = 6912;
N = n_particle;       % number of particles
n_conf = 8;

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

index_time = [2 3 6 11 12 15 20 21 24 29 30 33 38 39 42 47 48 51 56 57 60 65];
% index_time = 1:length(t);
% index_time = 21;

Temperature = [0.5:0.1:1.0 1.2 1.5 2.0 3.0 5.0];
% Temperature = [1.0];
%% start loop
for T_i = Temperature % loop over temperature
    tic
    
    iT = iT+1;
    disp(['T = ',num2str(T_i,'%.1f')])
    
    fraction_c = zeros(length(index_time),1);
    for ic = 1:n_conf % loop over configuration
        disp(ic)
        %% time origin
        %%%%%%%%%%%%%%%% load rg.mat file %%%%%%%%%%%%%%%%
        filename_rg = ['../rg_mat/rg_py_',num2str(T_i,'%.1f'),'.mat'];
        load(filename_rg)
        
        clear Vs_0
        it = 1;
        rg_all = rg_all_t_c(:,:,:,it,ic);
        [p_0,coord_0] = calculate_shape(rg_all);
        for i1 = 1:size(rg_all,3)
            [V,D] = eig(rg_all(:,:,i1));
            [lambda,ind] = sort(diag(D));
            Vs_0(:,:,i1) = V(:,ind);
        end
        
        %%%%%%%%%%%%%%%% load trajectory file %%%%%%%%%%%%%%%%
        filename = ['../',num2str(T_i,'%.1f'),'/eq.',num2str(ic,'%.0f'),'.dump'];
        dump = readmatrix(filename,'FileType','text');
        
        index = find(isnan(dump(:,3))==0);
        
        id = dump(index,1);
        r = dump(index,3:5);
        l = dump(1:3,1:2);
        
        type = dump(index,2);
        
        %% PBC
        Lx = (l(1,2)-l(1,1));
        Ly = (l(2,2)-l(2,1));
        Lz = (l(3,2)-l(3,1));
        L = [Lx, Ly, Lz];     % box size
        pbc = [1, 1, 1];      % boundary conditions
        r = r-round(r./L).*L;
        
        itime = 0;
        %% loop over time
        for it = index_time
            disp(['t = ', num2str(t(it))])
            itime = itime+1;
            rg_all = rg_all_t_c(:,:,:,it,ic);
            [p_t,coord_t] = calculate_shape(rg_all);
            for i1 = 1:size(rg_all,3)
                [V,D] = eig(rg_all(:,:,i1));
                [lambda,ind] = sort(diag(D));
                Vs(:,:,i1) = V(:,ind);
            end
            
            %%%%%%%%%%%%%%%% calculate correlation %%%%%%%%%%%%%%%%
            Ct = (3*(dot(Vs_0(:,3,:),Vs(:,3,:))).^2-1)/2;
            Ct = Ct(:);
            C_th = ((1-(1/2/pi))^2*3-1)/2;
%             C_th = 0.6;
            Ct_center = Ct-C_th;
            fraction(itime,1) = sum(Ct_center>0)/N;
            t_MC(itime,1) = t(it);

        end
        fraction_c = fraction_c + fraction;
    end
    fraction_c = fraction_c/n_conf;
    
    fraction_T(:,iT) = fraction_c;
       
    toc
end


%% Figures
close all
% linecolor

color_parula = flipud(parula(100));
index_color = round(0.44./Temperature*100);

color_order = color_parula(index_color,:);

% load('YlOrRd.mat')
% color = rgba(:,1:3);
% index_color = round((0.44./Temperature*0.75+0.25)*length(color));
% color_order = color(index_color,:);

figure(2)
load('fraction.mat')
hold on
box on
set(gca, 'ColorOrder', color_order)
t_MC = t(index_time)
plot(t_MC/1000,fraction_T,'LineWidth',2)

load('tau_M.mat')
tau_M = tau_M(2:end)/1000;
for i = 1:length(Temperature)
    F_tau = griddedInterpolant(t_MC/1000,fraction_T(:,i),'linear');
    fraction_tau_M(i) = F_tau(tau_M(i));
    plot(tau_M(i),fraction_tau_M(i),'o','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor',color_order(i,:),...
        'MarkerSize',8)
end

load susceptibility_1.2_1.0.mat
for i = 1:length(Temperature)
    F_max = griddedInterpolant(t_MC/1000,fraction_T(:,i),'linear');
    fraction_tmax(i) = F_max(tmax(i)/1000);
    plot(tmax(i)/1000,fraction_tmax(i),'^','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor',color_order(i,:),...
        'MarkerSize',8)
end

% plot([1 1e9]/1000,[0 0],'-','Color','#666666')

set(gca, 'XScale', 'log')
xlim([1e-3 1e4])
ylim([0 1])
xticks(10.^[-3:1:4])

xlabel('\it{t}','FontSize',28)

set(gca,'LineWidth',2)
set(gcf,'Position',[200,100,600,600])
set(gca,'position',[0.18    0.22   0.7376    0.7376])
set(gca,'FontSize',28,'FontName','Arial')

ylabel('Fraction','FontSize',26)

