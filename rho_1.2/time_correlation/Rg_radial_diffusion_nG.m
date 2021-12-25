clear;
close all;

n_conf = 4;
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
for T_i = fliplr([0.1 0.4 0.7 1.0 5.0]) % loop over temperature
    tic
    
    iT = iT+1;
    disp(['T = ',num2str(T_i,'%.1f')])
    
    % load rg.mat file
    filename_rg = ['../rg_mat/rg_py_',num2str(T_i,'%.1f'),'.mat'];
    load(filename_rg)
    
    % assign color representing the temperature
    cp = parula(20);
    
    %% preallocation
    MSDrjk = zeros(length(t),3);
    MQDrjk = zeros(length(t),3);
    MSDrjk_c = zeros(length(t),3);
    nG = zeros(length(t),3);
    nG_c = zeros(length(t),3);
    
    for ic = 1:n_conf % loop over configuration
        disp(ic)
        
        %% time origin
        clear Vs_0
        it = 1;
        rg_all = rg_all_t_c(:,:,:,it,ic);
        [p_0,coord_0] = calculate_shape(rg_all);
        %         for i1 = 1:size(rg_all,3)
        %             [V,D] = eig(rg_all(:,:,i1));
        %             [lambda,ind] = sort(diag(D));
        %             Vs_0(:,:,i1) = V(:,ind);
        %         end
        
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
        
        %%%%%%%%%%%%%%%% calculate pair distance and binning %%%%%%%%%%%%%%%%
%         rjk = zeros(N,N,3);
%         for i = 1:3
%             rjk(:,:,i) = r(:,i)-r(:,i)';
%             rjk(:,:,i) = rjk(:,:,i) - round(rjk(:,:,i) /L(i) ) * L_times_pbc(i);
%         end
%         djk = sqrt(sum(rjk.^2,3));
%         
%         rc = 1;
%         weight = exp(-djk.^2/(2*rc^2));       % weighted by gaussian distribution
%         weight_sum = sum(weight);                   % sum up the weight
%         
% %         index_bin = ceil(djk/dr);
%         clear djk        
%         
%         rjk = zeros(N,N,3);
%         for i = 1:3
%             rjk(:,:,i) = r(:,i)-r(:,i)';
% %             rjk(:,:,i) = rjk(:,:,i) - round(rjk(:,:,i) /L(i) ) * L_times_pbc(i);
%         end

        r0 = r;
        CM0 = mean(r);
        CMw0 = mean(rw);
        r0 = r0-CM0;
        %% loop over time
        for it = 2:n_frame
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
            
%             rjk_t = zeros(N,N,3);
%             for i = 1:3
%                 rjk_t(:,:,i) = r(:,i)-r(:,i)';
% %                 rjk_t(:,:,i) = rjk_t(:,:,i) - round(rjk_t(:,:,i) /L(i) ) * L_times_pbc(i);
%             end        
            
%             SDrjk = (sum((rjk_t-rjk).^2,3));
%             QDrjk = SDrjk.^2;
%             clear rjk_t

            rt = r;
            CMt = mean(r);
            CMwt = mean(rw);
            rt = rt-CMt;
            
            for is = 1:3
%                 index_shape = 1:length(state0);%==is;
                index_shape = state0==is;
                
%                 MSDrjk(it,is) = mean(sum(SDrjk(:,index_shape).*(weight(:,index_shape)./weight_sum(:,index_shape))));
%                 MQDrjk(it,is) = mean(sum(QDrjk(:,index_shape).*(weight(:,index_shape)./weight_sum(:,index_shape))));
                
                MSDrjk(it,is) = mean(sum((rt(index_shape,:)-r0(index_shape,:)).^2,2));
                MQDrjk(it,is) = mean((sum((rt(index_shape,:)-r0(index_shape,:)).^2,2)).^2);
                
                nG(it,is) = 3*MQDrjk(it,is)/(5*MSDrjk(it,is)^2)-1;
            end
            
            MSDrjk_c(it,:) =  MSDrjk_c(it,:) + MSDrjk(it,:);
            nG_c(it,:) = nG_c(it,:) + nG(it,:);
        end % end of time loop
               
    end % end of configuration loop
    MSDrjk_c = MSDrjk_c/n_conf;
    nG_c = nG_c/n_conf;
    
    toc
    

    figure(1)
    plot(t,nG_c(:,1:3),'-','LineWidth',2)
    set(gca, 'XScale', 'log')
%     set(gca, 'YScale', 'log')
    xlim([1e0 1e7])
%     ylim([0 1])
    box on
    
    xlabel('$t$','FontSize',24,'Interpreter','latex')
    ylabel('$\textrm{MSD}$','FontSize',24,'Interpreter','latex')
    set(gca,'LineWidth',2)
    set(gcf,'Position',[200,100,600,600])
    set(gca,'FontSize',28,'FontName','Times New Roman')
end% end of temperature loop