% Calculate tau_LC
clear;
close all;

n_conf = 1;
N = 6912;
dr = 0.25;
rc = 1.35;
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
% index_time = [2 3 6 11 12 15 20 21 24 29 30 33 38 39 42 47 48 51 56 57 60 65];
% index_time = 2;
index_time = 1:length(t);

%% preallocation
C_LC = zeros(length(t),6);
C_LC_c = zeros(length(t),6);
%% start loop
for T_i = fliplr([0.4:0.1:1.0 2.0 5.0]) % loop over temperature
    tic
    
    iT = iT+1;
    disp(['T = ',num2str(T_i,'%.1f')])
    
    % assign color to temperature
    cp = parula(50);
    color = cp(round(T_i*10),:);
       
    for ic = 1:n_conf % loop over configuration
        disp(ic)
        
        %% time origin
        it = 1;
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
        
        %%%%%%%%%%%%%%%% calculate pair distance and binning %%%%%%%%%%%%%%%%
        rjk = zeros(N,N,3);
        for i = 1:3
            rjk(:,:,i) = r(:,i)-r(:,i)';
            rjk(:,:,i) = rjk(:,:,i) - round(rjk(:,:,i) /L(i) ) * L_times_pbc(i);
        end
        
        djk = sqrt(sum(rjk.^2,3));
        
        bond_0 = zeros(N);
        bond_0(:) = logical(djk<rc);
        bond_0 = sparse(bond_0);
        
        clear djk
        clear rjk
        
        %% loop over time
        %         for it = 2:n_frame
        for it = index_time
            %%%%%%%%%%%%%%%% calculate pair distance and binning %%%%%%%%%%%%%%%%
            index_t = index((it-1)*N+(1:N));
            r = dump(index_t,3:5);
            
            rjk = zeros(N,N,3);
            for i = 1:3
                rjk(:,:,i) = r(:,i)-r(:,i)';
                rjk(:,:,i) = rjk(:,:,i) - round(rjk(:,:,i) /L(i) ) * L_times_pbc(i);
            end
            
            djk = sqrt(sum(rjk.^2,3));
            
            bond_t = zeros(N);
            bond_t(:) = logical(djk<rc);
            bond_t = sparse(bond_t);
            
            clear djk
            clear rjk
            
            %% compare connectivty
%             C_bond = 1-((bond_t.*bond_0)+((1-bond_t).*(1-bond_0))); % break and form
            C_bond = 0+(((1-bond_t).*bond_0)); % break
            C_LC(it,iT) = length(find(sum(C_bond)==0))/N;
            
        end % end of time loop
        C_LC_c(:,iT) = C_LC_c(:,iT) + C_LC(:,iT);
        
    end % end of configuration loop
    C_LC_c(:,iT) = C_LC_c(:,iT)/n_conf;
    
    [M,I_LC] = min(abs(C_LC_c-exp(-1)));
    t_LC = t(I_LC);
    
    toc
    
    %% figures
    figure(1)
    box on
    hold on
    plot(t(index_time),C_LC_c(index_time,iT),'LineWidth',2,'Color',color)
    
    % figure_settings
    set(gca, 'XScale', 'log')
    %     set(gca, 'YScale', 'log')
    xlim([1e0 1e7])
    xticks(10.^([1,4,7]))
    ylim([0 1])
    xlabel('$t$','FontSize',24,'Interpreter','latex')
    ylabel('$C_{LC}(t)$','FontSize',24,'Interpreter','latex')
    set(gca,'LineWidth',2)
    set(gcf,'Position',[200,100,600,600])
    set(gca,'FontSize',28,'FontName','Times New Roman')
    
end% end of temperature loop