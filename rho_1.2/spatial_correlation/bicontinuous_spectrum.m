clear
close all

n_particle = 6912;
N = n_particle;       % number of particles
n_conf = 8;
n_grid = 64;

%%
figure(1)
hold on
% % grids
% for iq = -5:3
%     plot([0.1 10],10^iq*([0.1 10]).^(-4),'-','Color','#666666')
% end
% 
% for iq = -5:3
%     plot([0.1 10],10^iq*([0.1 10]).^(-1),'-','Color','#666666')
% end

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
dr = 0.1;
rr = dr:dr:10;

%% temperature

Temperature = [0.5:0.1:1.0 1.2 1.5 2.0 3.0 5.0];
% Temperature = [1.0];

color = flipud(parula(100));
index_color = round(0.44./Temperature*100);

load tau_M.mat
tau_M = tau_M(2:end); % 0.5 to 5.0

%% start loop
for T_i = Temperature % loop over temperature
    tic
    
    iT = iT+1;
    disp(['T = ',num2str(T_i,'%.1f')])
    
    [min_dtau_M, index_time] = min(abs(t-tau_M(Temperature==T_i)));
    
    mean_MC_c = zeros(length(index_time),1);
    mean_GC_c = zeros(length(index_time),1);
    
    Pq_c = zeros(n_grid/2,1);
    
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
%         l = load('../l.txt');
        
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
            Ct_center = Ct-mean(Ct);
%             Ct_center = Ct-C_th;
            
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
            
            dl = Lx/n_grid; % pixel size
            [XX, YY, ZZ] = meshgrid((dl:dl:Lx)-Lx/2,(dl:dl:Ly)-Ly/2,(dl:dl:Lz)-Lz/2);
            F = scatteredInterpolant(Xs, Ys, Zs, ps, 'linear');
            CC = F(XX, YY, ZZ);
            
            SP = zeros(size(CC));
            
            [B,I_CC] = sort(CC(:));
            
            SP(I_CC(n_grid^3/2+1:end)) = 1;
%             SP(CC>0) = 1;
%             SP((XX.^2+YY.^2+ZZ.^2)<5) = 1;
            
            SP_FFT = fftn(SP,[n_grid,n_grid,n_grid]);
            
            abs_SP_FFT = abs(SP_FFT).^2;
            
            qx = 2*pi*(0:n_grid-1)/Lx;
            qy = 2*pi*(0:n_grid-1)/Ly;
            qz = 2*pi*(0:n_grid-1)/Lz;
            
            [qXX, qYY, qZZ] = meshgrid(qx,qy,qz);
            
            N_s = sum(SP,'all');
            V_s = Lx*Ly*Lz;
            
            %%
%             scatter3(qXX(:),qYY(:),qZZ(:),10,SP(:),'filled')
%             axis equal
%             
%             caxis([0 1])
            
            %%
            q = sqrt(qXX.^2+qYY.^2+qZZ.^2);
            
            dq = 2*pi/Lx;
            qq = dq*(1:n_grid/2);
            
            index_q = ceil(q/dq);
            
            Pq = zeros(length(qq),1);
            
            for iq = 1:n_grid/2
                vq = 4*pi*qq(iq)^2*dq/8;
                Pq(iq) = sum(abs_SP_FFT(index_q==iq))/vq/(N_s^2/2)*2;
                
            end
            
%             plot(qq,Pq/n_grid^3)
%             set(gca, 'XScale', 'log')
%             set(gca, 'YScale', 'log')

        end
        Pq_c = Pq_c + Pq;
        
    end
    Pq_c = Pq_c/n_conf;
    
    for i_r = 1:length(rr)
        Pr(i_r) = trapz(qq',qq'.*sin(qq'*rr(i_r)).*Pq_c)/rr(i_r);
    end
    
    %%
    figure(1)
    hold on
       
    plot(qq,Pq_c,'Color',color(index_color(iT),:),'LineWidth',2)
    set(gca, 'XScale', 'log')
    set(gca, 'YScale', 'log')
    
    xlabel('q','FontSize',24)
    ylabel('P(q)','FontSize',24)
    
    box on
    
    xlim([0.4 10])
    ylim([1e-4 1])
    
    set(gca,'LineWidth',2)
    set(gcf,'Position',[200,100,300,300])
    set(gca,'FontSize',24,'FontName','Arial')
    
    pause(0.1)
       
    toc
end


