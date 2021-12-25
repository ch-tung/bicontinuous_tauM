% Does those neighbor particlles around oblate clusters tend to stay around
% them?

clear;
close all;

n_conf = 1;
N = 6912;
dr = 0.25;
rmax = 5;
rr = (dr:dr:rmax)';
len_rr = length(rr);

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
index_time = [2 3 6 11 12 15 20 21 24 29 30 33 38 39 42 47 48 51 56 57 60 65];

%% figure background
figure(1)
hold on
% plot(t,t.^2'*(10.^(-20:2:20)),'-','Color','#C0C0C0')
% plot(t,t.^1'*(10.^(-20:2:20)),'-','Color','#C0C0C0')
% plot(t,t.^0.5'*(10.^(-20:2:20)),'-','Color','#C0C0C0')

%% start loop
for T_i = fliplr([5.0]) % loop over temperature
    tic
    
    iT = iT+1;
    disp(['T = ',num2str(T_i,'%.1f')])
    
    % load rg.mat file
    filename_rg = ['../rg_mat/rg_py_',num2str(T_i,'%.1f'),'.mat'];
    load(filename_rg)
    
    % assign color representing the temperature
    cp = parula(20);
    
    %% preallocation
    MDrjk_L = zeros(length(t),3);
    MSDrjk = zeros(length(t),3);
    MSDrjk_L = zeros(length(t),3);
    MSDrjk_T = zeros(length(t),3);
    
    MQDrjk = zeros(length(t),3);
    
    MDrjk_L_c = zeros(length(t),3);
    MSDrjk_c = zeros(length(t),3);
    MSDrjk_L_c = zeros(length(t),3);
    MSDrjk_T_c = zeros(length(t),3);
    
    nG = zeros(length(t),3);
    nG_c = zeros(length(t),3);
    
    for ic = 1:n_conf % loop over configuration
        disp(ic)
        
        %% time origin
        clear Vs_0
        it = 1;
        rg_all = rg_all_t_c(:,:,:,it,ic);
        [p_0,coord_0] = calculate_shape(rg_all);
        for i1 = 1:size(rg_all,3)
            [V,D] = eig(rg_all(:,:,i1));
            [lambda,ind] = sort(diag(D));
            Vs_0(:,:,i1) = V(:,ind);
        end
        u_0 = Vs_0(:,3,:);
        u_0 = permute(u_0,[2 3 1]);
        
        %%%%%%%%%%%%%%%% sort Rg into 3 type %%%%%%%%%%%%%%%%
        nstate = 3;
        state0 = zeros(1*N,1);
        [S_p,I_p] = sort(p_0(:,3));
        ni_s = floor(1*N/nstate);
        
        for is = 1
            state0(I_p(ni_s*(is-1)+1:ni_s*(is))) = is;
        end
        
        index_nonspere = find(state0==0);
        [S_S7,I_S7] = sort(coord_0(index_nonspere,2));
        for is = 2:3
            state0(index_nonspere(I_S7(ni_s*(is-2)+1:ni_s*(is-1)))) = is;
        end
        
        %%%%%%%%%%%%%%%% load trajectory file %%%%%%%%%%%%%%%%
        filename = ['../',num2str(T_i,'%.1f'),'/eq.',num2str(ic,'%.0f'),'.dump'];
        dump = readmatrix(filename,'FileType','text');
        
        index = find(isnan(dump(:,3))==0);
        index_t = index((it-1)*N+(1:N));
        
        id = dump(index_t,1);
        r = dump(index_t,3:5);
        l = dump(1:3,1:2);
                
        type = dump(index_t,2);
        
        % parameters from MD simulation
        Lx = (l(1,2)-l(1,1));
        Ly = (l(2,2)-l(2,1));
        Lz = (l(3,2)-l(3,1));
        L = [Lx, Ly, Lz];     % box size
        pbc = [1, 1, 1];      % boundary conditions
        L_times_pbc = L .* pbc; % deal with boundary conditions
        
        rw = r - floor((r-l(:,1)')./L).*L;
        r0 = r;
        CM0 = mean(r);
        CMw0 = mean(rw);
        r0 = r0-CM0;
        
%         r = r-CMw0;
        
        %%%%%%%%%%%%%%%% calculate pair distance and binning %%%%%%%%%%%%%%%%
        rjk = zeros(N,N,3);
        rjk_pbc = zeros(N,N,3);
        for i = 1:3
            rjk(:,:,i) = r(:,i)-r(:,i)';
            rjk_pbc(:,:,i) = rjk(:,:,i) - round(rjk(:,:,i) /L(i) ) * L_times_pbc(i);
%             rjk(:,:,i) = rjk(:,:,i) - round(rjk(:,:,i) /L(i) ) * L_times_pbc(i);
            rjk(:,:,i) = r(:,i)-0*r(:,i)';
        end
        
        djk = sqrt(sum(rjk_pbc.^2,3));
        
        ujk = rjk_pbc./djk;
        
        rc = 1;
        weight = exp(-djk.^2/(2*rc^2));       % weighted by gaussian distribution
        weight(weight<1e-8) = 0;
        weight = sparse(weight);
        weight_sum = sum(weight);             % sum up the weight
        
        weight_normalized = weight./weight_sum;
        
%         index_bin = ceil(djk/dr);
        clear djk     
        clear rjk_pbc 
        
        %% loop over time
%         for it = 2:n_frame
        for it = index_time
            rg_all = rg_all_t_c(:,:,:,it,ic);
            [p_t,coord_t] = calculate_shape(rg_all);
            %             for i1 = 1:size(rg_all,3)
            %                 [V,D] = eig(rg_all(:,:,i1));
            %                 [lambda,ind] = sort(diag(D));
            %                 Vs(:,:,i1) = V(:,ind);
            %             end
            %%%%%%%%%%%%%%%% calculate pair distance and binning %%%%%%%%%%%%%%%%
            index_t = index((it-1)*N+(1:N));
            r = dump(index_t,3:5);
            rw = r - floor((r-l(:,1)')./L).*L;
            rt = r;
            CMt = mean(r);
            CMwt = mean(rw);
            rt = rt-CMt;
            
%             r = r-CMwt;
            
            rjk_t = zeros(N,N,3);
            rjk_pbc_t = zeros(N,N,3);
            for i = 1:3
                rjk_t(:,:,i) = r(:,i)-0*r(:,i)';
%                 rjk_pbc_t(:,:,i) = rjk_t(:,:,i) - round(rjk_t(:,:,i) /L(i) ) * L_times_pbc(i);
                %             rjk(:,:,i) = rjk(:,:,i) - round(rjk(:,:,i) /L(i) ) * L_times_pbc(i);
            end
            
%             djk_t = sqrt(sum(rjk_pbc_t.^2,3));
            
%             ujk_t = rjk_pbc_t./djk_t;

            displacement = rjk_t-rjk;
            
            Drjk_L = (sum(displacement.*ujk,3));
            Drjk_L(boolean(eye(N))) = 0;  
            
            SDrjk = (sum(displacement.^2,3));
            SDrjk_L = (sum(displacement.*ujk,3).^2); % square of inner product
%             SDrjk_L = (sum((rjk_t-rjk).*u_0,3).^2);
            SDrjk_L(boolean(eye(N))) = 0;            
            SDrjk_T = SDrjk-SDrjk_L;
            SDrjk_T(boolean(eye(N))) = 0;  
            
            QDrjk = SDrjk.^2;
            clear rjk_t
            clear rjk_pbc_t
            clear djk_t
%             clear ujk
%             clear ujk_t            

            for is = 1:3
%                 index_shape = 1:length(state0);%==is;
                index_shape = state0==is;
                
                weight_index = weight_normalized(:,index_shape);
                
                MDrjk_L(it,is) = mean(sum(Drjk_L(:,index_shape).*weight_index));
                MSDrjk(it,is) = mean(sum(SDrjk(:,index_shape).*weight_index));
                MSDrjk_L(it,is) = mean(sum(SDrjk_L(:,index_shape).*weight_index));
                MSDrjk_T(it,is) = mean(sum(SDrjk_T(:,index_shape).*weight_index));
                
                MQDrjk(it,is) = mean(sum(QDrjk(:,index_shape).*weight_index));
                
%                 MSDrjk(it,is) = mean(mean(SDrjk(:,index_shape)));               
%                 MQDrjk(it,is) = mean(mean(QDrjk(:,index_shape)));
                
%                 MSDrjk(it,is) = mean(sum((rt(index_shape,:)-r0(index_shape,:)).^2,2));
%                 MQDrjk(it,is) = mean((sum((rt(index_shape,:)-r0(index_shape,:)).^2,2)).^2);
                
                nG(it,is) = 3*MQDrjk(it,is)/(5*MSDrjk(it,is)^2)-1;
            end
            
            MDrjk_L_c(it,:) =  MDrjk_L_c(it,:) + MDrjk_L(it,:);
            MSDrjk_c(it,:) =  MSDrjk_c(it,:) + MSDrjk(it,:);
            MSDrjk_L_c(it,:) =  MSDrjk_L_c(it,:) + MSDrjk_L(it,:);
            MSDrjk_T_c(it,:) =  MSDrjk_T_c(it,:) + MSDrjk_T(it,:);
            
            nG_c(it,:) = nG_c(it,:) + nG(it,:);
        end % end of time loop
               
    end % end of configuration loop
    MDrjk_L_c = MDrjk_L_c/n_conf;
    MSDrjk_c = MSDrjk_c/n_conf;
    MSDrjk_L_c = MSDrjk_L_c/n_conf;
    MSDrjk_T_c = MSDrjk_T_c/n_conf;
    
    nG_c = nG_c/n_conf;
    
    RMSDrjk_c = sqrt(MSDrjk_c);
    RMSDrjk_L_c = sqrt(MSDrjk_L_c);
    RMSDrjk_T_c = sqrt(MSDrjk_T_c);
    
    cos_theta = RMSDrjk_L_c./RMSDrjk_c;
    theta = acos(cos_theta);
    
    toc
    
    %%
    linestyle = {'--','-.',':'};
    newcolors_hex = {'0072BD','D95319','EDB120'};
    for i = 1:length(newcolors_hex)
        newcolors(i,:) = sscanf(newcolors_hex{i},'%2x%2x%2x',[1 3])/255;
    end
%     newcolors = parula(3);
    set(gca, 'ColorOrder', newcolors)

    figure(1)
    hold on
    for i_s = 1:3
        plot(t(index_time),nG_c(index_time,i_s),linestyle{i_s},'LineWidth',2)
    end
    set(gca, 'XScale', 'log')
%     set(gca, 'YScale', 'log')
    xlim([1e0 1e7])
%     ylim([1e0 1e7])
    box on
    
    xlabel('$t$','FontSize',24,'Interpreter','latex')
    ylabel('$\alpha_{2}$','FontSize',24,'Interpreter','latex')
    set(gca,'LineWidth',2)
    set(gcf,'Position',[200,100,600,600])
    set(gca,'FontSize',28,'FontName','Times New Roman')
    
    %%
    figure(2)
    set(gca, 'ColorOrder', newcolors)
    
    hold on
    for i_s = 1:3
        plot(t(index_time),MSDrjk_c(index_time,i_s),linestyle{i_s},'LineWidth',2)
        plot(t(index_time),MSDrjk_L_c(index_time,i_s),linestyle{i_s},'LineWidth',2)
        plot(t(index_time),MSDrjk_T_c(index_time,i_s),linestyle{i_s},'LineWidth',2)
    end
    
    set(gca, 'XScale', 'log')
    set(gca, 'YScale', 'log')
    xlim([1e0 1e7])
    ylim([1e-7 1e0])
    box on
    
    xlabel('$t$','FontSize',24,'Interpreter','latex')
    ylabel('$\textrm{MSD}$','FontSize',24,'Interpreter','latex')
    set(gca,'LineWidth',2)
    set(gcf,'Position',[200,100,600,600])
    set(gca,'FontSize',28,'FontName','Times New Roman')
    
    %%
    figure(3)
    hold on
    set(gca, 'ColorOrder', newcolors)
    
    for i_s = 1:3
        plot(t(index_time),theta((index_time),i_s),linestyle{i_s},'LineWidth',2)
    end
    
    set(gca, 'XScale', 'log')
%     set(gca, 'YScale', 'log')
    xlim([1e0 1e7])
    ylim([0 pi/2])
    yticks([0:pi/4:pi/2])
    yticklabels({'0','\pi/4','\pi/2'})
    
    box on
    
    xlabel('$t$','FontSize',24,'Interpreter','latex')
    ylabel('$\theta$','FontSize',24,'Interpreter','latex')
    set(gca,'LineWidth',2)
    set(gcf,'Position',[200,100,600,600])
    set(gca,'FontSize',28,'FontName','Times New Roman')
end% end of temperature loop