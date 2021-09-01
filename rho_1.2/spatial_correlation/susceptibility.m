clear;
close all;

n_conf = 8;
N = 6912;
dr = 0.25;
rmax = 5;
rr = (dr:dr:rmax)';
len_rr = length(rr);
a0 = 1;

load('tau_M.mat')
tau_M = tau_M(2:end);

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
%% figure background
figure(1)
hold on
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
% plot(t,t.^2'*(10.^(-20:2:20)),'-','Color','#C0C0C0')
% plot(t,t.^1'*(10.^(-20:2:20)),'-','Color','#C0C0C0')
% plot(t,t.^0.5'*(10.^(-20:2:20)),'-','Color','#C0C0C0')

%% start loop
T = [0.5:0.1:1.0 1.2 1.5 2.0 3.0 5.0];
% T = 1.0;

color_parula = flipud(parula(100));
index_color = round(0.44./T*100);
color_order = color_parula(index_color,:);

set(gca, 'ColorOrder', color_order)

chi_c_T = zeros(length(index_time),length(T));
for T_i = T % loop over temperature
    tic
    
    iT = iT+1;
    disp(['T = ',num2str(T_i,'%.1f')])
    
    % load rg.mat file
    filename_rg = ['../rg_mat/rg_py_',num2str(T_i,'%.1f'),'.mat'];
    load(filename_rg)
    
    % assign color representing the temperature
    cp = parula(20);
    
    %% preallocation

    
    for ic = 1:n_conf % loop over configuration
        disp(ic)
        
        %% time origin
        clear Vs_0
        it = 1;
        rg_all = rg_all_t_c(:,:,:,it,ic);
        
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
        
        %% loop over time
%         for it = 2:n_frame
        for it = index_time
            index_t = index((it-1)*N+(1:N));
            r = dump(index_t,3:5);
            rw = r - floor((r-l(:,1)')./L).*L;
            rt = r;
            CMt = mean(r);
            CMwt = mean(rw);
            rt = rt-CMt;
            
            Drt = sqrt(sum((r0-rt).^2,2));
            Qt = zeros(size(Drt));
            Qt(Drt<=a0) = 1;
            chi(it) = var(Qt);
        end % end of time loop
        chi_c_T(:,iT) = chi_c_T(:,iT) + chi';       
        
    end % end of configuration loop
    chi_c_T(:,iT) = chi_c_T(:,iT)/n_conf;
    chi_c_T(:,iT) = chi_c_T(:,iT);
       
    [max_chi_c_T, index_tmax] = max(chi_c_T(:,iT));
    tmax(iT,1) = t(index_time(index_tmax));
    
    toc
    
    plot(t(index_time),chi_c_T(:,iT),'-','LineWidth',2)
end% end of temperature loop

%%
t_fine = 1.1.^(0:169);
chi_c_T_interp = zeros(length(t_fine),length(T));
iT_p = 0;
for T_i = T % loop over temperature
    iT_p = iT_p+1;

    F = griddedInterpolant(log(t(index_time(2:length(index_time))))',chi_c_T(2:end,iT_p),'spline');
    chi_c_T_interp(:,iT_p) = F(log(t_fine));
    
    plot(t_fine,chi_c_T_interp(:,iT_p),'-','LineWidth',2)
    
    [max_chi_c_T_i, index_tmax_i] = max(chi_c_T_interp(:,iT_p));
    tmax(iT_p,1) = t_fine(index_tmax_i);
    plot(tmax(iT_p,1),chi_c_T_interp(index_tmax_i,iT_p),'o')
end
