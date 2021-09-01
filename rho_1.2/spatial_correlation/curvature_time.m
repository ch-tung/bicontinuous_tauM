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

% index_time = [2 3 6 11 12 15 20 21 24 29 30 33 38 39 42 47 48 51 56 57 60 65];
index_time = 1:length(t);
% index_time = 21;

Temperature = [0.5:0.1:1.0 1.2 1.5 2.0 3.0 5.0];
% Temperature = [1.0];
%% start loop
for T_i = Temperature % loop over temperature
    tic
    
    iT = iT+1;
    disp(['T = ',num2str(T_i,'%.1f')])
    
    mean_MC_c = zeros(length(index_time),1);
    mean_GC_c = zeros(length(index_time),1);
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
%             Ct_center = Ct-mean(Ct);
            Ct_center = Ct-C_th;
            
            p = [Ct_center];
            
            %         % scatter
            %         figure(1)
            %         scatter3(X,Y,Z,50,p(index_selection),'filled')
            %         axis equal
            
            %% contour
            index_t0 = ((1-1)*N+(1:N));
            index_t = ((it-1)*N+(1:N));
            rt0 = r(index_t0,:);
            rt = r(index_t,:);
            
            %%%%%%%%%%%%%%%% define region %%%%%%%%%%%%%%%%
            thickness = 16;
            center = 0;
            
            origin_X = 0;
            origin_Y = 0;
            origin_Z = 0;
            
            %         index_selection = find(abs(rt0(:,3)-center)<thickness/2);
            index_selection = 1:N;
            X = rt0(index_selection,1);
            Y = rt0(index_selection,2);
            Z = rt0(index_selection,3);
            
            % shadow
            Xs = [];
            Ys = [];
            Zs = [];
            ps = [];
            for ix = -1:1
                for iy = -1:1
                    for iz = -1:1
                        Xs = [Xs;X+Lx*ix];
                        Ys = [Ys;Y+Ly*iy];
                        Zs = [Zs;Z+Lz*iz];
                        ps = [ps;p];
                    end
                end
            end
            
            edge = 2.5;
            index_shadow = find(abs((Xs)+origin_X)<abs(l(1,2)-l(1,1))/2+edge&...
                abs((Ys)+origin_Y)<abs(l(2,2)-l(2,1))/2+edge);
            
            Xs = Xs(index_shadow);
            Ys = Ys(index_shadow);
            Zs = Zs(index_shadow);
            ps = ps(index_shadow,:);
            
            dl = Lx/10; % pixel size
            n_edge = 1;
            [XX, YY, ZZ] = meshgrid((-n_edge*dl:dl:Lx+n_edge*dl)-Lx/2,(-n_edge*dl:dl:Ly+n_edge*dl)-Ly/2,(-n_edge*dl:dl:Lz+n_edge*dl)-Ly/2);
            F = scatteredInterpolant(Xs, Ys, Zs, ps, 'natural');
            CC = F(XX, YY, ZZ);
            
            %% isosurface
            [f, v] = isosurface(XX, YY, ZZ, CC, 0);
            
            if length(v)<100
                mean_GC(itime,1) = NaN;
                mean_MC(itime,1) = NaN;
                t_MC(itime,1) = t(it);
                continue
            end
            
            [GC, MC, Area] = curvatures(v(:,1),v(:,2),v(:,3),f);
            
            i_interior = sum(min(abs(v),[Lx/2,Lx/2,Lx/2])==[Lx/2,Lx/2,Lx/2],2)==0;
            
%             MC_i = MC(i_interior);
%             MC_i = MC_i(abs(MC_i)<5);
%             mean_MC(itime,1) = mean(MC_i);
            
            MC_i = MC(i_interior);
            A_i = Area(i_interior);
            mean_MC(itime,1) = sum(MC_i.*A_i')/sum(A_i);
            
%             GC_i = GC(i_interior);
%             GC_i = GC_i(abs(GC_i)<5);
%             mean_GC(itime,1) = mean(GC_i);

            GC_i = GC(i_interior);
            mean_GC(itime,1) = sum(GC_i.*A_i')/sum(A_i);
            
            t_MC(itime,1) = t(it);
            
%             figure(1)
%             histogram(MC_i,-5:0.2:5,'Normalization','pdf')
        end
        mean_MC_c = mean_MC_c + mean_MC;
        mean_GC_c = mean_GC_c + mean_GC;
    end
    mean_MC_c = mean_MC_c/n_conf;
    mean_GC_c = mean_GC_c/n_conf;
    
    mean_MC_T(:,iT) = mean_MC_c;
    mean_GC_T(:,iT) = mean_GC_c;
       
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
hold on
box on
set(gca, 'ColorOrder', color_order)

plot(t_MC/1000,mean_MC_T,'LineWidth',2)

load('tau_M.mat')
tau_M = tau_M(2:end)/1000;
for i = 1:length(Temperature)
    F_tau = griddedInterpolant(t_MC/1000,mean_MC_T(:,i),'linear');
    MC_tau_M(i) = F_tau(tau_M(i));
    plot(tau_M(i),MC_tau_M(i),'o','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor',color_order(i,:),...
        'MarkerSize',12)
end

plot([1 1e9]/1000,[0 0],'-','Color','#666666')

set(gca, 'XScale', 'log')
xlim([10 1e5]/1000)
ylim([-0.5 0.5])

xlabel('$t$','FontSize',28,'Interpreter','latex')
ylabel('$\left\langle H \right\rangle$','FontSize',28,'Interpreter','latex')
set(gca,'LineWidth',2)
set(gcf,'Position',[200,100,800,600])
set(gca,'FontSize',28,'FontName','Times New Roman')

%%
figure(3)
hold on
box on
set(gca, 'ColorOrder', color_order)

mean_GC_T(mean_GC_T>0) = 0;

plot(t_MC,2*pi./sqrt(-6*mean_GC_T),'LineWidth',2)
set(gca, 'XScale', 'log')

for i = 1:length(Temperature)
    F_tau = griddedInterpolant(t_MC,mean_GC_T(:,i),'linear');
    GC_tau_M(i) = F_tau(tau_M(i));
    plot(tau_M(i),2*pi./sqrt(-6*GC_tau_M(i)),'o','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor',color_order(i,:),...
        'MarkerSize',12)
end

xlim([10 1e6])
ylim([0 20])

xlabel('$t$','FontSize',24,'Interpreter','latex')
ylabel('$2\pi\left\langle k^2 \right\rangle^{-\frac{1}{2}}$','FontSize',24,'Interpreter','latex')
set(gca,'LineWidth',2)
set(gcf,'Position',[200,100,800,600])
set(gca,'FontSize',28,'FontName','Times New Roman')
